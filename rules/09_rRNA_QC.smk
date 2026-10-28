# Report-only rRNA content QC: bowtie2 of (subsampled) trimmed reads vs rRNA FASTA.
# Default FASTA is extracted from the run GTF (rRNA / Mt_rRNA / …); override with rrna_fasta.


rule extract_rrna_fasta:
    """Pull rRNA(/Mt_rRNA) transcripts from GTF+genome when rrna_fasta is not set."""
    input:
        genome = genome,
        gtf = gtf,
    output:
        fasta = os.path.join(rrna_ref_dir, "rrna_from_gtf.fa"),
        gtf_filt = os.path.join(rrna_ref_dir, "rrna_from_gtf.gtf"),
    params:
        biotypes = ",".join(rrna_biotypes),
        filter_script = os.path.join(workflow.basedir, "scripts", "filter_rrna_gtf.py"),
    log:
        os.path.join(work_dir, log_dir, "extract_rrna_fasta.log"),
    conda:
        "../envs/rrna.yml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
        nodes=1,
        cpus_per_task=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fasta})" "$(dirname {log})"
        echo "[LOG] Filtering GTF for biotypes: {params.biotypes}" | tee {log}

        # External script (not an indented heredoc) — avoids Python IndentationError
        # when Snakemake expands the shell block.
        python3 {params.filter_script} \
            --gtf {input.gtf} \
            --output {output.gtf_filt} \
            --biotypes "{params.biotypes}" \
            >> {log} 2>&1

        gffread -w {output.fasta} -g {input.genome} {output.gtf_filt} >> {log} 2>&1
        if [ ! -s {output.fasta} ]; then
            echo "ERROR: extracted rRNA FASTA is empty after gffread." | tee -a {log}
            echo "Filtered GTF line count: $(wc -l < {output.gtf_filt})" | tee -a {log}
            echo "Set config rrna_fasta to an external rRNA FASTA instead." | tee -a {log}
            exit 1
        fi
        echo "[LOG] Extracted $(grep -c '^>' {output.fasta} || true) rRNA sequences" | tee -a {log}
        """


# Index input: custom FASTA from config, otherwise GTF-extracted file.
_rrna_fa_for_index = (
    rrna_fasta
    if rrna_fasta not in (None, False, "")
    else rules.extract_rrna_fasta.output.fasta
)


rule bowtie2_index_rrna:
    input:
        fasta = _rrna_fa_for_index,
    output:
        flag = touch(os.path.join(rrna_index_dir, "index.done")),
    params:
        prefix = os.path.join(rrna_index_dir, "rrna"),
    log:
        os.path.join(work_dir, log_dir, "bowtie2_index_rrna.log"),
    conda:
        "../envs/rrna.yml"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
        nodes=1,
        cpus_per_task=4,
    shell:
        r"""
        set -euo pipefail
        mkdir -p {rrna_index_dir} "$(dirname {log})"
        echo "[LOG] Building bowtie2 rRNA index from {input.fasta}" | tee {log}
        bowtie2-build --threads {threads} {input.fasta} {params.prefix} >> {log} 2>&1
        if [ ! -e {params.prefix}.1.bt2 ] && [ ! -e {params.prefix}.1.bt2l ]; then
            echo "ERROR: bowtie2 index missing after build" | tee -a {log}
            exit 1
        fi
        echo "✓ rRNA bowtie2 index ready at {params.prefix}" | tee -a {log}
        """


rule rrna_bowtie2:
    input:
        read1 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R1_trimmed.fastq.gz"),
        read2 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R2_trimmed.fastq.gz"),
        status = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.status"),
        index_flag = os.path.join(rrna_index_dir, "index.done"),
    output:
        bt2_log = os.path.join(work_dir, rrna_dir, "{sample}", "{sample}.bowtie2.log"),
        summary = os.path.join(work_dir, rrna_dir, "{sample}", "{sample}.rrna_summary.txt"),
    params:
        prefix = os.path.join(rrna_index_dir, "rrna"),
        subsample = rrna_subsample_reads,
        outdir = os.path.join(work_dir, rrna_dir, "{sample}"),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
        nodes=1,
        cpus_per_task=4,
    log:
        os.path.join(work_dir, log_dir, "rrna_bowtie2_{sample}.log"),
    conda:
        "../envs/rrna.yml"
    shell:
        r"""
        set -euo pipefail

        if [ "$(head -n 1 {input.status})" != "ok" ]; then
            echo "ERROR: refusing rRNA QC for skipped sample {wildcards.sample}" | tee -a {log}
            exit 1
        fi

        mkdir -p {params.outdir} "$(dirname {log})"
        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT

        echo "[LOG] Subsampling {params.subsample} pairs for {wildcards.sample}" | tee {log}
        seqtk sample -s 42 {input.read1} {params.subsample} > "$tmpdir/R1.fastq"
        seqtk sample -s 42 {input.read2} {params.subsample} > "$tmpdir/R2.fastq"

        if [ ! -s "$tmpdir/R1.fastq" ] || [ ! -s "$tmpdir/R2.fastq" ]; then
            echo "ERROR: subsampled FASTQs empty for {wildcards.sample}" | tee -a {log}
            exit 1
        fi

        echo "[LOG] bowtie2 vs rRNA index {params.prefix}" | tee -a {log}
        bowtie2 \
            -x {params.prefix} \
            -1 "$tmpdir/R1.fastq" \
            -2 "$tmpdir/R2.fastq" \
            -p {threads} \
            --no-unal \
            -S /dev/null \
            2> {output.bt2_log}

        if ! grep -qi 'overall alignment rate' {output.bt2_log}; then
            echo "ERROR: bowtie2 log missing overall alignment rate" | tee -a {log}
            cat {output.bt2_log} >> {log} || true
            exit 1
        fi

        rate=$(grep -i 'overall alignment rate' {output.bt2_log} | tail -n1 | awk '{{print $1}}' | tr -d '%')
        echo -e "sample\trrna_alignment_rate_pct\n{wildcards.sample}\t$rate" > {output.summary}
        echo "✓ rRNA bowtie2 {wildcards.sample}: ${{rate}}% aligned" | tee -a {log}
        cat {output.bt2_log} >> {log}
        """


rule rrna_qc:
    input:
        logs = lambda wildcards: expand(
            os.path.join(work_dir, rrna_dir, "{sample}", "{sample}.bowtie2.log"),
            sample=get_kept_samples(wildcards),
        ),
        salmon_meta = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "aux_info", "meta_info.json"),
            sample=get_kept_samples(wildcards),
        ),
        kept = os.path.join(work_dir, "metrics", "kept_samples_after_trim.txt"),
    output:
        tsv = os.path.join(work_dir, "metrics", "rrna_qc_mqc.tsv"),
    params:
        rrna_root = os.path.join(work_dir, rrna_dir),
        salmon_root = os.path.join(work_dir, salmon_dir),
        script = os.path.join(workflow.basedir, "scripts", "rrna_qc.py"),
    log:
        os.path.join(work_dir, log_dir, "rrna_qc.log"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
        nodes=1,
        cpus_per_task=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.tsv})"
        echo "[LOG] Starting rrna_qc $(date -Is)" > {log}

        manifest="$(dirname {output.tsv})/rrna_qc_manifest.tsv"
        : > "$manifest"
        while IFS= read -r sample || [ -n "$sample" ]; do
            sample=$(echo "$sample" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            [ -z "$sample" ] && continue
            btlog="{params.rrna_root}/$sample/$sample.bowtie2.log"
            meta="{params.salmon_root}/$sample/aux_info/meta_info.json"
            printf '%s\t%s\t%s\n' "$sample" "$btlog" "$meta" >> "$manifest"
        done < {input.kept}

        echo "[LOG] Manifest lines: $(wc -l < "$manifest")" >> {log}
        python3 {params.script} \
            --manifest "$manifest" \
            --output {output.tsv} \
            >> {log} 2>&1

        echo "[LOG] Done rrna_qc $(date -Is)" >> {log}
        """
