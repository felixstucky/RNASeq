#!/usr/bin/env python3
"""Aggregate bowtie2 rRNA alignment rates for MultiQC (report-only)."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Optional


COLUMNS = [
    "sample",
    "salmon_percent_mapped",
    "rrna_alignment_rate_pct",
    "rrna_reads_aligned",
    "rrna_reads_total",
    "bowtie2_log",
]

# bowtie2 stderr: "55.00% overall alignment rate"
OVERALL_RE = re.compile(
    r"([0-9.]+)%\s+overall alignment rate",
    re.IGNORECASE,
)
# PE: "X (Y%) aligned concordantly ..." plus discordant; overall rate is the headline metric.
ALIGNED_CONCORDANT_RE = re.compile(
    r"(\d+)\s+\([^)]+\)\s+aligned concordantly\s+exactly 1 time",
    re.IGNORECASE,
)
ALIGNED_DISCORDANT_RE = re.compile(
    r"(\d+)\s+\([^)]+\)\s+aligned discordantly 1 time",
    re.IGNORECASE,
)
PAIRS_RE = re.compile(r"^(\d+)\s+reads; of these:", re.IGNORECASE | re.MULTILINE)


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


def parse_bowtie2_log(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="replace")
    overall = OVERALL_RE.search(text)
    if not overall:
        raise ValueError(f"could not parse overall alignment rate in {path}")

    pairs = PAIRS_RE.search(text)
    total = int(pairs.group(1)) if pairs else ""

    # Best-effort aligned count (concordant exactly 1 + discordant); optional.
    aligned = ""
    if pairs:
        conc = ALIGNED_CONCORDANT_RE.findall(text)
        disc = ALIGNED_DISCORDANT_RE.search(text)
        # Also count ">1 time"
        conc_multi = re.search(
            r"(\d+)\s+\([^)]+\)\s+aligned concordantly\s+>1 times",
            text,
            re.IGNORECASE,
        )
        n = 0
        for m in conc:
            n += int(m)
        if conc_multi:
            n += int(conc_multi.group(1))
        if disc:
            n += int(disc.group(1))
        # PE log counts pairs; report pairs aligned
        aligned = n

    return {
        "rrna_alignment_rate_pct": float(overall.group(1)),
        "rrna_reads_aligned": aligned,
        "rrna_reads_total": total,
    }


def write_multiqc_tsv(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header_lines = [
        "# id: rrna_qc",
        "# section_name: 'rRNA QC (bowtie2)'",
        "# description: 'Fraction of (subsampled) trimmed reads aligning to rRNA references. Report-only — explains Salmon unmapped / high-GC libraries.'",
        "# format: 'tsv'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#   id: 'rrna_qc_table'",
        "#   title: 'rRNA QC'",
        "#   namespace: ''",
        "\t".join(COLUMNS),
    ]
    with output_path.open("w") as handle:
        handle.write("\n".join(header_lines) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in COLUMNS) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV: sample, bowtie2.log, [salmon meta_info.json]",
    )
    parser.add_argument("--output", required=True, help="MultiQC custom-content TSV")
    args = parser.parse_args()

    rows: list[dict[str, Any]] = []
    with Path(args.manifest).open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"ERROR: bad manifest line: {line}", file=sys.stderr)
                return 1
            sample, log_path = parts[0], Path(parts[1])
            meta_path = Path(parts[2]) if len(parts) > 2 and parts[2] else None
            if not log_path.is_file():
                print(f"ERROR: missing bowtie2 log for {sample}: {log_path}", file=sys.stderr)
                return 1
            parsed = parse_bowtie2_log(log_path)
            rows.append(
                {
                    "sample": sample,
                    "salmon_percent_mapped": load_salmon_percent_mapped(meta_path),
                    "rrna_alignment_rate_pct": parsed["rrna_alignment_rate_pct"],
                    "rrna_reads_aligned": parsed["rrna_reads_aligned"],
                    "rrna_reads_total": parsed["rrna_reads_total"],
                    "bowtie2_log": str(log_path),
                }
            )

    if not rows:
        print("ERROR: no samples in rRNA QC manifest", file=sys.stderr)
        return 1

    write_multiqc_tsv(rows, Path(args.output))
    rates = [float(r["rrna_alignment_rate_pct"]) for r in rows]
    print(
        f"Wrote rRNA QC table: {args.output} "
        f"(n={len(rows)}, rrna% min={min(rates):.1f} median={sorted(rates)[len(rates)//2]:.1f} max={max(rates):.1f})",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
