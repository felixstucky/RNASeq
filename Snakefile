include: "rules/00_Variable_Setup.smk"
include: "rules/01_Symlink_Files.smk"
include: "rules/02_FastQC.smk"
include: "rules/03_Trimming.smk"
include: "rules/04_MultiQC.smk"
include: "rules/05_Salmon_Quant.smk"
include: "rules/07_Salmon_Mapping_QC.smk"
include: "rules/06_Generate_Count_Matrix.smk"
include: "rules/08_Kraken.smk"
include: "rules/09_rRNA_QC.smk"

# Tiny / orchestration jobs: run on the submit host, not via SLURM srun.
# Avoids cpus-per-task mismatches and conda/singularity hangs on QC aggregates.
localrules: all, create_links, select_trimmed_samples, salmon_mapping_qc, kraken_species_qc, rrna_qc

SAMPLES = sample
wildcard_constraints:
    sample = "[^/]+"

rule all:
    input:
        transcript_counts = os.path.join(
            work_dir, count_matrix_dir, "TranscriptCountMatrix.csv"
        ),
        transcript_tpm = os.path.join(
            work_dir, count_matrix_dir, "TranscriptTPMMatrix.csv"
        ),
        count_Matrix = os.path.join(
            work_dir, count_matrix_dir, "GeneCountMatrix.csv"
        ),
        count_Matrix_non_zero = os.path.join(
            work_dir, count_matrix_dir, "GeneCountMatrix_non_zero_count.csv"
        ),
        count_Matrix_protein_coding_non_zero = os.path.join(
            work_dir, count_matrix_dir, "GeneCountMatrix_non_zero_count_protein_coding_only.csv"
        ),
        salmon_quant = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "quant.sf"),
            sample=get_kept_samples(wildcards),
        ),
        fastqc_html = expand(
            os.path.join(work_dir, "metrics", "{sample}_R{read}_fastqc.html"),
            sample=SAMPLES,
            read=["1", "2"],
        ),
        fastqc_zip = expand(
            os.path.join(work_dir, "metrics", "{sample}_R{read}_fastqc.zip"),
            sample=SAMPLES,
            read=["1", "2"],
        ),
        trimmed_fastqc_html = lambda wildcards: expand(
            os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.html"),
            sample=get_kept_samples(wildcards),
            read=["1", "2"],
        ),
        trimmed_fastqc_zip = lambda wildcards: expand(
            os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.zip"),
            sample=get_kept_samples(wildcards),
            read=["1", "2"],
        ),
        trimmed_fastq = expand(
            os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R{read}_trimmed.fastq.gz"),
            sample=SAMPLES,
            read=["1", "2"],
        ),
        bbduk_stats = expand(
            os.path.join(work_dir, trimmed_fastq_dir, "{sample}_bbduk_stats.txt"),
            sample=SAMPLES,
        ),
        kept_after_trim = os.path.join(work_dir, "metrics", "kept_samples_after_trim.txt"),
        skipped_after_trim = os.path.join(work_dir, "metrics", "skipped_after_trim.tsv"),
        salmon_mapping_qc = os.path.join(
            work_dir, "metrics", "salmon_mapping_qc_mqc.tsv"
        ),
        kept_after_mapping = os.path.join(
            work_dir, "metrics", "kept_samples_after_mapping.txt"
        ),
        skipped_after_mapping = os.path.join(
            work_dir, "metrics", "skipped_after_mapping.tsv"
        ),
        kraken_reports = lambda wildcards: expand(
            os.path.join(work_dir, kraken_dir, "{sample}", "{sample}.bracken.kreport2"),
            sample=get_kept_samples(wildcards),
        ),
        kraken_species_qc = os.path.join(
            work_dir, "metrics", "kraken_species_qc_mqc.tsv"
        ),
        rrna_qc = os.path.join(work_dir, "metrics", "rrna_qc_mqc.tsv"),
        multiqc_html = os.path.join(work_dir, "multiqc", "multiqc_report.html")
