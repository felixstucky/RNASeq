#!/usr/bin/env python3
"""Aggregate Kraken2 reports for contamination QC (report-only by default)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional


TAXID_UNCLASSIFIED = 0
TAXID_BACTERIA = 2
TAXID_ARCHAEA = 2157
TAXID_VIRUSES = 10239
TAXID_HUMAN = 9606

COLUMNS = [
    "sample",
    "salmon_percent_mapped",
    "unclassified_pct",
    "bacteria_pct",
    "viral_pct",
    "archaea_pct",
    "human_pct",
    "expected_taxid",
    "expected_pct",
    "top_taxon",
    "top_taxid",
    "top_pct",
]


def parse_threshold(raw: str) -> Optional[float]:
    value = str(raw).strip().lower()
    if value in {"", "null", "none", "false", "na", "nan"}:
        return None
    return float(raw)


def parse_kreport(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            rows.append(
                {
                    "pct": float(parts[0]),
                    "clade_reads": int(float(parts[1])),
                    "tax_reads": int(float(parts[2])),
                    "rank": parts[3].strip(),
                    "taxid": int(parts[4]),
                    "name": parts[5].strip(),
                }
            )
    return rows


def clade_pct(rows: list[dict[str, Any]], taxid: int) -> float:
    for row in rows:
        if row["taxid"] == taxid:
            return row["pct"]
    return 0.0


def load_salmon_percent_mapped(path: Optional[Path]) -> Any:
    if path is None or not path.is_file():
        return ""
    with path.open() as handle:
        meta = json.load(handle)
    percent_mapped = meta.get("percent_mapped")
    if percent_mapped is None and meta.get("num_processed"):
        num_mapped = float(meta.get("num_mapped", 0))
        num_processed = float(meta["num_processed"])
        percent_mapped = (100.0 * num_mapped / num_processed) if num_processed else 0.0
    return percent_mapped if percent_mapped is not None else ""


def summarize_report(
    sample: str,
    report_path: Path,
    expected_taxid: int,
    salmon_meta_path: Optional[Path] = None,
) -> dict[str, Any]:
    rows = parse_kreport(report_path)
    if not rows:
        raise ValueError(f"empty or unreadable kreport: {report_path}")

    unclassified_pct = clade_pct(rows, TAXID_UNCLASSIFIED)
    for row in rows:
        if row["rank"] == "U" or row["name"].lower() == "unclassified":
            unclassified_pct = row["pct"]
            break

    top_row = None
    for row in rows:
        if row["taxid"] in {TAXID_UNCLASSIFIED} or row["rank"] == "U":
            continue
        if row["rank"].startswith("S"):
            if top_row is None or row["pct"] > top_row["pct"]:
                top_row = row
    if top_row is None:
        for row in rows:
            if row["taxid"] in {TAXID_UNCLASSIFIED} or row["rank"] == "U":
                continue
            if top_row is None or row["pct"] > top_row["pct"]:
                top_row = row

    return {
        "sample": sample,
        "salmon_percent_mapped": load_salmon_percent_mapped(salmon_meta_path),
        "unclassified_pct": unclassified_pct,
        "bacteria_pct": clade_pct(rows, TAXID_BACTERIA),
        "viral_pct": clade_pct(rows, TAXID_VIRUSES),
        "archaea_pct": clade_pct(rows, TAXID_ARCHAEA),
        "human_pct": clade_pct(rows, TAXID_HUMAN),
        "expected_taxid": expected_taxid,
        "expected_pct": clade_pct(rows, expected_taxid),
        "top_taxon": top_row["name"] if top_row else "",
        "top_taxid": top_row["taxid"] if top_row else "",
        "top_pct": top_row["pct"] if top_row else "",
    }


def write_multiqc_tsv(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header_lines = [
        "# id: kraken_species_qc",
        "# section_name: 'Kraken Contamination QC'",
        "# description: 'Kraken2 contamination profile (unclassified / bacteria / viral) vs Salmon mapping'",
        "# format: 'tsv'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#   id: 'kraken_species_qc_table'",
        "#   title: 'Kraken Contamination QC'",
        "#   namespace: ''",
        "\t".join(COLUMNS),
    ]
    with output_path.open("w") as handle:
        handle.write("\n".join(header_lines) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in COLUMNS) + "\n")


def check_expected_pct(
    rows: list[dict[str, Any]],
    min_expected_pct: Optional[float],
) -> list[str]:
    if min_expected_pct is None:
        return []
    failures = []
    for row in rows:
        value = row.get("expected_pct")
        if value == "" or value is None:
            failures.append(f"{row['sample']}: missing expected_pct")
            continue
        if float(value) < float(min_expected_pct):
            failures.append(
                f"{row['sample']}: expected_pct={value} < min_expected_pct={min_expected_pct}"
                f" (bacteria={row.get('bacteria_pct')}%, top={row.get('top_taxon')})"
            )
    return failures


def load_manifest(path: Path) -> list[tuple[str, Path, Optional[Path]]]:
    rows = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"manifest line {line_no}: need sample\\treport[\\tsalmon_meta]")
            sample = parts[0].strip()
            report = Path(parts[1].strip())
            meta = Path(parts[2].strip()) if len(parts) > 2 and parts[2].strip() else None
            rows.append((sample, report, meta))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV with columns: sample, kreport2, [salmon meta_info.json]",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-taxid", type=int, default=10090)
    parser.add_argument("--min-expected-pct", default="null")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(f"[LOG] Reading manifest: {args.manifest}", flush=True)
    entries = load_manifest(Path(args.manifest))
    print(f"[LOG] {len(entries)} sample(s) in manifest", flush=True)

    rows = []
    for sample, report_path, meta_path in entries:
        if not report_path.is_file():
            print(f"ERROR: missing report for {sample}: {report_path}", file=sys.stderr)
            return 1
        rows.append(
            summarize_report(sample, report_path, args.expected_taxid, meta_path)
        )

    output_path = Path(args.output)
    write_multiqc_tsv(rows, output_path)
    print(f"Wrote Kraken contamination QC table: {output_path}", flush=True)

    threshold = parse_threshold(args.min_expected_pct)
    failures = check_expected_pct(rows, threshold)
    if failures:
        print("ERROR: Kraken expected-taxon QC failed:", file=sys.stderr)
        for line in failures:
            print(f"  - {line}", file=sys.stderr)
        return 1

    if threshold is None:
        print("Report-only mode: expected-taxon fail gate disabled", flush=True)
    else:
        print(
            f"All samples passed min_expected_pct={threshold} for taxid {args.expected_taxid}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
