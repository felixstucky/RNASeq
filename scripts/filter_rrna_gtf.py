#!/usr/bin/env python3
"""Filter a GTF for rRNA-related biotypes (used by extract_rrna_fasta)."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


BIOTYPE_RE = re.compile(r'(?:gene_biotype|transcript_biotype)\s+"([^"]+)"')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--biotypes",
        required=True,
        help="Comma-separated biotype names to keep",
    )
    args = parser.parse_args()

    biotypes = {x.strip() for x in args.biotypes.split(",") if x.strip()}
    if not biotypes:
        print("ERROR: empty --biotypes", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    n_lines = 0
    n_exon = 0
    n_transcript = 0
    n_gene = 0
    with args.gtf.open() as hin, args.output.open("w") as hout:
        for line in hin:
            if line.startswith("#"):
                continue
            found = BIOTYPE_RE.findall(line)
            if not any(b in biotypes for b in found):
                continue
            hout.write(line)
            n_lines += 1
            feature = line.split("\t")[2] if "\t" in line else ""
            if feature == "exon":
                n_exon += 1
            elif feature == "transcript":
                n_transcript += 1
            elif feature == "gene":
                n_gene += 1

    print(
        f"[LOG] Wrote {n_lines} GTF feature lines "
        f"(gene={n_gene}, transcript={n_transcript}, exon={n_exon}) -> {args.output}",
        flush=True,
    )
    if n_lines == 0:
        print(
            "ERROR: no GTF features matched rrna_biotypes="
            + str(sorted(biotypes))
            + ". Set config rrna_fasta to an external rRNA FASTA instead.",
            file=sys.stderr,
        )
        return 1
    if n_exon == 0:
        print(
            "ERROR: matched rRNA biotypes but found 0 exon lines. "
            "gffread needs exons/CDSs to write a FASTA. "
            "Set config rrna_fasta to an external rRNA FASTA "
            "(e.g. SILVA / Ensembl ncRNA) instead.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
