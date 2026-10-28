rule kraken2_classify:
    input:
        read1 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R1_trimmed.fastq.gz"),
        read2 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R2_trimmed.fastq.gz"),
        status = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.status"),
    output:
        report = os.path.join(work_dir, kraken_dir, "{sample}", "{sample}.kreport2"),
        kraken = os.path.join(work_dir, kraken_dir, "{sample}", "{sample}.kraken2"),
        bracken = os.path.join(work_dir, kraken_dir, "{sample}", "{sample}.bracken"),
        bracken_report = os.path.join(
            work_dir, kraken_dir, "{sample}", "{sample}.bracken.kreport2"
        ),
    params:
        db = kraken_db,
        subsample = kraken_subsample_reads,
        read_length = read_length,
        outdir = os.path.join(work_dir, kraken_dir, "{sample}"),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
        nodes=1,
        cpus_per_task=8,
    log:
        os.path.join(work_dir, log_dir, "kraken2_{sample}.log")
    conda:
        "../envs/kraken.yml"
    shell:
        r"""
        set -euo pipefail

        if [ "$(head -n 1 {input.status})" != "ok" ]; then
            echo "ERROR: refusing Kraken for skipped sample {wildcards.sample}" | tee -a {log}
            exit 1
        fi

        mkdir -p {params.outdir}
        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT

        echo "[LOG] Subsampling {params.subsample} pairs for {wildcards.sample}" | tee -a {log}
        seqtk sample -s 42 {input.read1} {params.subsample} > "$tmpdir/R1.fastq"
        seqtk sample -s 42 {input.read2} {params.subsample} > "$tmpdir/R2.fastq"

        if [ ! -s "$tmpdir/R1.fastq" ] || [ ! -s "$tmpdir/R2.fastq" ]; then
            echo "ERROR: subsampled FASTQs are empty for {wildcards.sample}" | tee -a {log}
            exit 1
        fi

        echo "[LOG] Running Kraken2 for {wildcards.sample}" | tee -a {log}
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --paired \
            --report {output.report} \
            --output {output.kraken} \
            "$tmpdir/R1.fastq" \
            "$tmpdir/R2.fastq" \
            >> {log} 2>&1

        if [ ! -s {output.report} ]; then
            echo "ERROR: Kraken2 report missing/empty for {wildcards.sample}" | tee -a {log}
            exit 1
        fi

        echo "[LOG] Running Bracken (species) for {wildcards.sample}" | tee -a {log}
        bracken \
            -d {params.db} \
            -i {output.report} \
            -o {output.bracken} \
            -w {output.bracken_report} \
            -r {params.read_length} \
            -l S \
            >> {log} 2>&1

        if [ ! -s {output.bracken} ] || [ ! -s {output.bracken_report} ]; then
            echo "ERROR: Bracken outputs missing/empty for {wildcards.sample}" | tee -a {log}
            exit 1
        fi

        echo "✓ Kraken2/Bracken complete for {wildcards.sample}" | tee -a {log}
        """


rule kraken_species_qc:
    input:
        reports = lambda wildcards: expand(
            os.path.join(work_dir, kraken_dir, "{sample}", "{sample}.kreport2"),
            sample=get_kept_samples(wildcards),
        ),
        salmon_meta = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "aux_info", "meta_info.json"),
            sample=get_kept_samples(wildcards),
        ),
        kept = os.path.join(work_dir, "metrics", "kept_samples_after_trim.txt"),
    output:
        tsv = os.path.join(work_dir, "metrics", "kraken_species_qc_mqc.tsv"),
    params:
        expected_taxid = expected_taxid,
        min_expected_pct = (
            "null" if min_expected_pct in (None, False) else min_expected_pct
        ),
        kraken_root = os.path.join(work_dir, kraken_dir),
        salmon_root = os.path.join(work_dir, salmon_dir),
        script = os.path.join(workflow.basedir, "scripts", "kraken_species_qc.py"),
    log:
        os.path.join(work_dir, log_dir, "kraken_species_qc.log"),
    # stdlib-only script: avoid conda/singularity startup hangs on tiny QC jobs
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1,
        cpus_per_task=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.tsv})"
        echo "[LOG] Starting kraken_species_qc $(date -Is)" > {log}

        manifest="$(dirname {output.tsv})/kraken_species_qc_manifest.tsv"
        : > "$manifest"
        while IFS= read -r sample || [ -n "$sample" ]; do
            sample=$(echo "$sample" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            [ -z "$sample" ] && continue
            report="{params.kraken_root}/$sample/$sample.kreport2"
            meta="{params.salmon_root}/$sample/aux_info/meta_info.json"
            printf '%s\t%s\t%s\n' "$sample" "$report" "$meta" >> "$manifest"
        done < {input.kept}

        echo "[LOG] Manifest lines: $(wc -l < "$manifest")" >> {log}
        python3 {params.script} \
            --manifest "$manifest" \
            --output {output.tsv} \
            --expected-taxid {params.expected_taxid} \
            --min-expected-pct {params.min_expected_pct} \
            >> {log} 2>&1

        echo "[LOG] Done kraken_species_qc $(date -Is)" >> {log}
        """
