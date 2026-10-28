#!/usr/bin/env python3
"""Aggregate Salmon mapping metrics; exclude low-mapping samples from count matrices."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Optional


COLUMNS = [
    "sample",
    "mapping_qc_status",
    "percent_mapped",
    "min_mapping_rate",
    "num_processed",
    "num_mapped",
    "num_decoy_fragments",
    "num_compatible_targets",
    "library_types",
    "expected_format",
]


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


_LIB_META_KEYS = {
    "read_files",
    "expected_format",
    "compatible_fragment_ratio",
    "num_compatible_fragments",
    "num_assigned_fragments",
    "num_frags_with_concordant_consistent_mappings",
    "num_frags_with_inconsistent_or_orphan_mappings",
    "strand_mapping_bias",
}


def _format_library_types(lib_data: Optional[dict[str, Any]]) -> str:
    if not lib_data:
        return ""
    format_counts = {
        key: value
        for key, value in lib_data.items()
        if key not in _LIB_META_KEYS and isinstance(value, (int, float))
    }
    if format_counts:
        ranked = sorted(format_counts.items(), key=lambda item: item[1], reverse=True)
        return ";".join(f"{name}:{count}" for name, count in ranked[:5])
    expected = lib_data.get("expected_format")
    return str(expected) if expected is not None else ""


def parse_threshold(raw: str) -> Optional[float]:
    value = str(raw).strip().lower()
    if value in {"", "null", "none", "false", "na", "nan"}:
        return None
    return float(raw)


def parse_sample(
    sample: str,
    meta_path: Path,
    lib_path: Optional[Path],
    min_mapping_rate: Optional[float],
) -> tuple[dict[str, Any], bool, str]:
    meta = _load_json(meta_path)
    lib_data = _load_json(lib_path) if lib_path and lib_path.is_file() else None

    percent_mapped = meta.get("percent_mapped")
    if percent_mapped is None and meta.get("num_processed"):
        num_mapped = float(meta.get("num_mapped", 0))
        num_processed = float(meta["num_processed"])
        percent_mapped = (100.0 * num_mapped / num_processed) if num_processed else 0.0

    passed = True
    reason = "pass"
    if min_mapping_rate is None:
        status = "pass"
    elif percent_mapped is None or percent_mapped == "":
        passed = False
        status = "fail"
        reason = "missing percent_mapped"
    elif float(percent_mapped) < float(min_mapping_rate):
        passed = False
        status = "fail"
        reason = f"percent_mapped={percent_mapped} < min_mapping_rate={min_mapping_rate}"
    else:
        status = "pass"

    row = {
        "sample": sample,
        "mapping_qc_status": status,
        "percent_mapped": percent_mapped if percent_mapped is not None else "",
        "min_mapping_rate": min_mapping_rate if min_mapping_rate is not None else "",
        "num_processed": meta.get("num_processed", ""),
        "num_mapped": meta.get("num_mapped", ""),
        "num_decoy_fragments": meta.get("num_decoy_fragments", ""),
        "num_compatible_targets": meta.get(
            "num_compatible_targets",
            meta.get("num_targets", ""),
        ),
        "library_types": _format_library_types(lib_data),
        "expected_format": (lib_data or {}).get("expected_format", ""),
    }
    return row, passed, reason


def write_multiqc_tsv(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header_lines = [
        "# id: salmon_mapping_qc",
        "# section_name: 'Salmon Mapping QC'",
        "# description: 'Per-sample Salmon mapping metrics (fail = excluded from count matrices)'",
        "# format: 'tsv'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#   id: 'salmon_mapping_qc_table'",
        "#   title: 'Salmon Mapping QC'",
        "#   namespace: ''",
        "\t".join(COLUMNS),
    ]
    with output_path.open("w") as handle:
        handle.write("\n".join(header_lines) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in COLUMNS) + "\n")


def write_skipped_mqc(skipped: list[tuple[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header_lines = [
        "# id: salmon_mapping_excluded",
        "# section_name: 'Excluded After Mapping QC'",
        "# description: 'Samples excluded from count matrices due to low Salmon mapping rate'",
        "# format: 'tsv'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#   id: 'salmon_mapping_excluded_table'",
        "#   title: 'Excluded After Mapping QC'",
        "#   namespace: ''",
        "sample\treason",
    ]
    with output_path.open("w") as handle:
        handle.write("\n".join(header_lines) + "\n")
        for sample, reason in skipped:
            handle.write(f"{sample}\t{reason}\n")


def load_manifest(path: Path) -> list[tuple[str, Path, Optional[Path]]]:
    rows = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"manifest line {line_no}: need sample\\tmeta[\\tlib]")
            sample = parts[0].strip()
            meta = Path(parts[1].strip())
            lib = Path(parts[2].strip()) if len(parts) > 2 and parts[2].strip() else None
            rows.append((sample, meta, lib))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV with columns: sample, meta_info.json, [lib_format_counts.json]",
    )
    parser.add_argument("--output", required=True, help="MultiQC mapping QC table")
    parser.add_argument("--kept-output", required=True)
    parser.add_argument("--skipped-output", required=True)
    parser.add_argument("--skipped-mqc-output", default="")
    parser.add_argument("--min-mapping-rate", default="30")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(f"[LOG] Reading manifest: {args.manifest}", flush=True)
    entries = load_manifest(Path(args.manifest))
    print(f"[LOG] {len(entries)} sample(s) in manifest", flush=True)

    threshold = parse_threshold(args.min_mapping_rate)
    rows = []
    kept = []
    skipped: list[tuple[str, str]] = []

    for sample, meta_path, lib_path in entries:
        if not meta_path.is_file():
            print(f"ERROR: missing meta_info.json for {sample}: {meta_path}", file=sys.stderr)
            return 1
        row, passed, reason = parse_sample(sample, meta_path, lib_path, threshold)
        rows.append(row)
        if passed:
            kept.append(sample)
        else:
            skipped.append((sample, reason))

    write_multiqc_tsv(rows, Path(args.output))

    kept_path = Path(args.kept_output)
    kept_path.parent.mkdir(parents=True, exist_ok=True)
    with kept_path.open("w") as handle:
        handle.write("\n".join(kept) + ("\n" if kept else ""))

    skipped_path = Path(args.skipped_output)
    with skipped_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "reason"])
        writer.writerows(skipped)

    if args.skipped_mqc_output:
        write_skipped_mqc(skipped, Path(args.skipped_mqc_output))

    print(f"Wrote Salmon mapping QC table: {args.output}", flush=True)
    print(f"Kept {len(kept)}/{len(entries)} samples for count matrices -> {kept_path}", flush=True)
    if skipped:
        print(f"Excluded {len(skipped)} sample(s) after mapping QC -> {skipped_path}", flush=True)
        for sample, reason in skipped:
            print(f"  - {sample}: {reason}", flush=True)

    if not kept:
        print("ERROR: No samples left after mapping QC filtering.", file=sys.stderr)
        return 1

    if threshold is None:
        print("min_mapping_rate disabled: all samples with Salmon quant kept", flush=True)
    else:
        print(f"Mapping QC threshold={threshold}% (exclude-only; pipeline continues)", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
