import csv
import os
import sys

gtf= config["gtf"]
genome= config["genome"]
genome_version = config["genome_version"]
metadata_file = config["metadata_file"]
read_length=config["read_length"]
sample = config["sample"]

work_dir = config["work_dir"]
genome_index_dir = config["genome_index_dir"]
salmon_kmer = config.get("salmon_kmer", 31)
min_mapping_rate = config.get("min_mapping_rate", 30)
min_fastq_bytes = int(config.get("min_fastq_bytes", 1000))
kraken_db = config.get(
    "kraken_db",
    "/data/cephfs-2/unmirrored/projects/cubit/current/static_data/app_support/kraken/2/minikraken_8GB_20200312",
)
kraken_subsample_reads = int(config.get("kraken_subsample_reads", 1000000))
# Mouse taxid; note Minikraken often lacks mouse and only has human + microbes.
expected_taxid = int(config.get("expected_taxid", 10090))
# Report-only by default (null). Accept legacy min_human_pct key if set.
min_expected_pct = config.get("min_expected_pct", config.get("min_human_pct", None))

# rRNA QC (bowtie2 vs rRNA FASTA). Report-only; explains Salmon unmapped fraction.
# If rrna_fasta is unset, refs are extracted from GTF biotypes (rRNA / Mt_rRNA / etc.).
rrna_fasta = config.get("rrna_fasta", None)
rrna_subsample_reads = int(config.get("rrna_subsample_reads", 1000000))
rrna_biotypes = config.get(
    "rrna_biotypes",
    ["rRNA", "Mt_rRNA", "rRNA_pseudogene", "ribonuclease_MRP_RNA", "ribozyme"],
)
if isinstance(rrna_biotypes, str):
    rrna_biotypes = [x.strip() for x in rrna_biotypes.split(",") if x.strip()]

for file_path, file_name in [
    (gtf, "GTF"),
    (genome, "Genome FASTA"),
    (metadata_file, "Metadata"),
    (kraken_db, "Kraken2 database"),
]:
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"ERROR: {file_name} file not found at {file_path}")

if rrna_fasta not in (None, False, "") and not os.path.exists(rrna_fasta):
    raise FileNotFoundError(f"ERROR: rRNA FASTA not found at {rrna_fasta}")

if not os.path.isdir(work_dir):
    os.makedirs(work_dir, exist_ok=True)
    print(f"✓ Created work directory: {work_dir}")

two_pass = config["two_pass"]
sequencing_strand = config["sequencing_strand"]


def normalize_sample_name(value):
    sample_value = str(value)
    if sample_value.endswith("quant.sf"):
        return os.path.basename(os.path.dirname(sample_value))
    return os.path.basename(os.path.normpath(sample_value))


def normalize_read_label(value):
    read_value = str(value).strip().upper()
    if read_value in {"R1", "1", "READ1", "READ_1"}:
        return "R1"
    if read_value in {"R2", "2", "READ2", "READ_2"}:
        return "R2"
    return read_value


def load_metadata_paths(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        required_cols = {"Sample_Name", "Read", "Absolute_Path"}
        if not reader.fieldnames or not required_cols.issubset(set(reader.fieldnames)):
            missing = sorted(required_cols - set(reader.fieldnames or []))
            raise ValueError(f"ERROR: Metadata file missing required columns: {missing}")

        paths = {}
        for row in reader:
            sample_name = normalize_sample_name(row["Sample_Name"])
            read_label = normalize_read_label(row["Read"])
            paths.setdefault(sample_name, {})[read_label] = row["Absolute_Path"]
        return paths


def filter_samples_by_fastq_size(samples, metadata_paths, min_bytes):
    kept = []
    skipped = []
    for sample_name in samples:
        reads = metadata_paths.get(sample_name, {})
        r1 = reads.get("R1")
        r2 = reads.get("R2")
        reasons = []

        if not r1 or not r2:
            missing = [label for label, path in (("R1", r1), ("R2", r2)) if not path]
            reasons.append(f"missing metadata paths for {','.join(missing)}")
        else:
            for label, path in (("R1", r1), ("R2", r2)):
                if not os.path.exists(path):
                    reasons.append(f"{label} not found: {path}")
                else:
                    size = os.path.getsize(path)
                    if size < min_bytes:
                        reasons.append(f"{label} too small ({size} < {min_bytes} bytes): {path}")

        if reasons:
            skipped.append((sample_name, "; ".join(reasons)))
        else:
            kept.append(sample_name)
    return kept, skipped


raw_sample = config["sample"]
if isinstance(raw_sample, (list, tuple)):
    all_samples = [normalize_sample_name(item) for item in raw_sample]
else:
    all_samples = [normalize_sample_name(raw_sample)]

metadata_paths = load_metadata_paths(metadata_file)
sample, skipped_samples = filter_samples_by_fastq_size(
    all_samples, metadata_paths, min_fastq_bytes
)

metrics_dir = os.path.join(work_dir, "metrics")
os.makedirs(metrics_dir, exist_ok=True)
skipped_samples_tsv = os.path.join(metrics_dir, "skipped_samples.tsv")
with open(skipped_samples_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample", "reason"])
    for sample_name, reason in skipped_samples:
        writer.writerow([sample_name, reason])

if skipped_samples:
    print(
        f"⚠ Skipping {len(skipped_samples)} sample(s) with missing/too-small FASTQs "
        f"(min_fastq_bytes={min_fastq_bytes}). See {skipped_samples_tsv}"
    )
    for sample_name, reason in skipped_samples:
        print(f"  - {sample_name}: {reason}")

if not sample:
    raise ValueError(
        "ERROR: No samples left after FASTQ size filtering. "
        f"Check metadata paths and min_fastq_bytes={min_fastq_bytes}."
    )

print(f"✓ Keeping {len(sample)}/{len(all_samples)} sample(s) for analysis")


log_dir = "Log"
symlink_dir = "Symlinks"
reference_dir = "Reference"
salmon_dir = "Salmon"
salmon_ref_dir = os.path.join(work_dir, reference_dir, f"Salmon_{genome_version}_k{salmon_kmer}")
salmon_index_dir = os.path.join(genome_index_dir, f"Salmon_{genome_version}_k{salmon_kmer}")
count_matrix_dir = "Count_Matrix"
trimmed_fastq_dir = "trimmed_fastqs"
kraken_dir = "Kraken"
rrna_dir = "rRNA"
rrna_ref_dir = os.path.join(work_dir, reference_dir, "rRNA")
rrna_index_dir = os.path.join(genome_index_dir, f"rRNA_{genome_version}")
# Effective FASTA used for indexing/alignment (custom override or extracted).
rrna_fasta_effective = (
    rrna_fasta
    if rrna_fasta not in (None, False, "")
    else os.path.join(rrna_ref_dir, "rrna_from_gtf.fa")
)
