rule multiqc:
    input:
        raw_fastqc = lambda wildcards: expand(
            os.path.join(work_dir, "metrics", "{sample}_R{read}_fastqc.zip"),
            sample=sample,
            read=["1", "2"],
        ),
        trimmed_fastqc = lambda wildcards: expand(
            os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.zip"),
            sample=get_kept_samples(wildcards),
            read=["1", "2"],
        ),
        bbduk_stats = lambda wildcards: expand(
            os.path.join(work_dir, trimmed_fastq_dir, "{sample}_bbduk_stats.txt"),
            sample=get_kept_samples(wildcards),
        ),
        salmon_quant = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "quant.sf"),
            sample=get_kept_samples(wildcards),
        ),
        meta_info = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "aux_info", "meta_info.json"),
            sample=get_kept_samples(wildcards),
        ),
        salmon_mapping_qc = os.path.join(work_dir, "metrics", "salmon_mapping_qc_mqc.tsv"),
        skipped_mapping_mqc = os.path.join(
            work_dir, "metrics", "salmon_mapping_excluded_mqc.tsv"
        ),
        skipped_after_mapping = os.path.join(
            work_dir, "metrics", "skipped_after_mapping.tsv"
        ),
        kraken_reports = lambda wildcards: expand(
            os.path.join(work_dir, kraken_dir, "{sample}", "{sample}.kreport2"),
            sample=get_kept_samples(wildcards),
        ),
        bracken_reports = lambda wildcards: expand(
            os.path.join(work_dir, kraken_dir, "{sample}", "{sample}.bracken.kreport2"),
            sample=get_kept_samples(wildcards),
        ),
        kraken_species_qc = os.path.join(work_dir, "metrics", "kraken_species_qc_mqc.tsv"),
        rrna_qc = os.path.join(work_dir, "metrics", "rrna_qc_mqc.tsv"),
    output:
        html = os.path.join(work_dir, "multiqc", "multiqc_report.html"),
        data = directory(os.path.join(work_dir, "multiqc", "multiqc_report_data"))
    params:
        search_dirs = (
            f"{work_dir}/metrics {work_dir}/metrics/trimmed "
            f"{work_dir}/{trimmed_fastq_dir} {work_dir}/{salmon_dir} "
            f"{work_dir}/{kraken_dir} {work_dir}/{rrna_dir}"
        ),
        outdir = os.path.join(work_dir, "multiqc"),
        title = "RNA-seq QC Report"
    log:
        os.path.join("logs", "multiqc", "multiqc.log")
    conda:
        "../envs/multiqc.yml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1,
        cpus_per_task=1,
    shell:
        """
        mkdir -p {params.outdir}

        multiqc \
            {params.search_dirs} \
            -o {params.outdir} \
            --filename multiqc_report.html \
            --title "{params.title}" \
            --force \
            --verbose \
            &> {log}

        echo "Syncing filesystem..." >> {log}
        sync
        sleep 15

        echo "Checking for output files..." >> {log}
        if [ ! -f {output.html} ]; then
            echo "ERROR: HTML report not found after sync" >> {log}
            echo "Directory contents:" >> {log}
            ls -lh {params.outdir}/ >> {log}
            exit 1
        fi

        if [ ! -d {output.data} ]; then
            echo "ERROR: Data directory not found after sync" >> {log}
            echo "Directory contents:" >> {log}
            ls -lh {params.outdir}/ >> {log}
            exit 1
        fi

        if [ ! -s {output.html} ]; then
            echo "ERROR: HTML report is empty" >> {log}
            exit 1
        fi

        echo "✓ MultiQC outputs verified successfully" >> {log}
        exit 0
        """
