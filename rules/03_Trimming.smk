rule trim_reads:
    input:
        r1 = os.path.join(work_dir, symlink_dir, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(work_dir, symlink_dir, "{sample}_R2.fastq.gz")
    output:
        r1 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R1_trimmed.fastq.gz"),
        r2 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R2_trimmed.fastq.gz"),
        stats = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_bbduk_stats.txt"),
        status = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.status"),
        flag = touch(os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.done"))
    params:
        ref = "adapters",
        extra = "ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen=75",
        # Leave headroom under the SLURM cgroup limit for the JVM itself.
        java_mem_mb = lambda wildcards, resources: max(1000, int(resources.mem_mb * 0.75)),
    threads: 2
    resources:
        # BBDuk + gzip I/O needs headroom beyond -Xmx; 8G was often OOM-killed on SLURM.
        mem_mb=16000,
        runtime=240,
        nodes=1,
        # Must match threads so srun --cpus-per-task does not exceed the SLURM allocation
        # when --default-resources includes cpus_per_task=1.
        cpus_per_task=2,
    log:
        os.path.join("logs", "trim_reads", "{sample}.log")
    conda:
        "../envs/bbduk.yml"
    shell:
        r"""
        skip_sample() {{
            reason="$1"
            echo "SKIP: {wildcards.sample}: $reason" | tee -a {log}
            gzip -c /dev/null > {output.r1}
            gzip -c /dev/null > {output.r2}
            printf "skipped\t%s\n" "$reason" > {output.stats}
            echo "skip" > {output.status}
            touch {output.flag}
            exit 0
        }}

        fail_hard() {{
            reason="$1"
            echo "ERROR: {wildcards.sample}: $reason" | tee -a {log}
            echo "---- last 80 log lines ----" >> {log}
            tail -n 80 {log} >> {log} 2>/dev/null || true
            exit 1
        }}

        mkdir -p "$(dirname {log})" "$(dirname {output.r1})"
        echo "[LOG] Starting trim_reads for {wildcards.sample}" > {log}

        # Lightweight input checks only (full gunzip -t on large Ceph FASTQs
        # can exhaust walltime before bbduk even starts).
        for fq in {input.r1} {input.r2}; do
            if [ ! -s "$fq" ]; then
                skip_sample "input missing or empty: $fq"
            fi
            # Confirm gzip magic bytes (1f 8b); do not scan the whole file.
            magic=$(od -An -tx1 -N2 "$fq" 2>> {log} | tr -d ' \n')
            if [ "$magic" != "1f8b" ]; then
                skip_sample "input is not gzip-compressed (bad magic): $fq"
            fi
        done

        echo "[LOG] Running bbduk with -Xmx{params.java_mem_mb}m (job mem_mb={resources.mem_mb})" >> {log}

        # Write to node-local scratch first. Concurrent gzip writers on Ceph often
        # die with "Stream closed" after processing completes (output may look fine
        # in the summary but BBDuk still exits non-zero).
        scratch_root="${{SLURM_TMPDIR:-${{TMPDIR:-/tmp}}}}"
        tmpdir="$scratch_root/bbduk_{wildcards.sample}_$$"
        mkdir -p "$tmpdir"
        trap 'rm -rf "$tmpdir"' EXIT

        tmp_r1="$tmpdir/{wildcards.sample}_R1_trimmed.fastq.gz"
        tmp_r2="$tmpdir/{wildcards.sample}_R2_trimmed.fastq.gz"
        tmp_stats="$tmpdir/{wildcards.sample}_bbduk_stats.txt"

        echo "[LOG] BBDuk temp dir: $tmpdir" >> {log}
        bbduk.sh \
            -Xmx{params.java_mem_mb}m \
            in1={input.r1} \
            in2={input.r2} \
            out1="$tmp_r1" \
            out2="$tmp_r2" \
            ref={params.ref} \
            {params.extra} \
            threads={threads} \
            stats="$tmp_stats" \
            overwrite=t \
            >> {log} 2>&1
        bbduk_status=$?

        result_reads=$(awk '/^Result:/{{print $2}}' {log} | tail -n 1)
        echo "[LOG] BBDuk exit=$bbduk_status Result reads=${{result_reads:-NA}}" | tee -a {log}

        if [ "$bbduk_status" -ne 0 ]; then
            # Clear FASTQ corruption signals only — do NOT treat generic Java
            # "Stream closed" / Exception as a bad library (common on shared FS).
            if grep -qiE 'unexpected end of (file|zlib)|Corrupt input|Invalid gzip|Premature EOF|Error reading' {log}; then
                skip_sample "bbduk failed with exit code $bbduk_status (likely bad FASTQ; see log)"
            fi
            # Salvage: summary shows retained reads and temp outputs exist.
            if [ -n "$result_reads" ] && [ "$result_reads" -gt 0 ] && [ -s "$tmp_r1" ] && [ -s "$tmp_r2" ] && [ -s "$tmp_stats" ]; then
                echo "WARN: bbduk exit $bbduk_status after non-zero Result; salvaging local outputs" | tee -a {log}
            else
                fail_hard "bbduk failed with exit code $bbduk_status (tool/infra error, not skipped)"
            fi
        fi

        echo "Validating trimmed outputs..." | tee -a {log}

        if [ ! -s "$tmp_r1" ] || [ ! -s "$tmp_r2" ]; then
            fail_hard "trimmed FASTQ missing/empty in scratch after bbduk"
        fi

        if [ ! -s "$tmp_stats" ]; then
            fail_hard "bbduk stats empty or missing in scratch after bbduk"
        fi

        if [ -z "$result_reads" ]; then
            fail_hard "could not parse BBDuk Result read count from log"
        fi
        if [ "$result_reads" -eq 0 ]; then
            skip_sample "bbduk succeeded but retained zero reads (empty after trim/minlen)"
        fi

        echo "[LOG] Copying outputs from scratch to {work_dir}" >> {log}
        cp -f "$tmp_r1" {output.r1}
        cp -f "$tmp_r2" {output.r2}
        cp -f "$tmp_stats" {output.stats}

        echo "✓ Trimming and validation complete for {wildcards.sample}" | tee -a {log}
        echo "ok" > {output.status}
        touch {output.flag}
        exit 0
        """


checkpoint select_trimmed_samples:
    input:
        status = expand(
            os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.status"),
            sample=sample,
        ),
    output:
        kept = os.path.join(work_dir, "metrics", "kept_samples_after_trim.txt"),
        skipped = os.path.join(work_dir, "metrics", "skipped_after_trim.tsv"),
    log:
        os.path.join(work_dir, log_dir, "select_trimmed_samples.log"),
    run:
        kept = []
        skipped_rows = []
        for sample_name, status_path in zip(sample, input.status):
            with open(status_path) as handle:
                status_value = handle.read().strip().splitlines()[0].strip().lower()
            if status_value == "ok":
                kept.append(sample_name)
            else:
                stats_path = os.path.join(
                    work_dir, trimmed_fastq_dir, f"{sample_name}_bbduk_stats.txt"
                )
                reason = "trim skipped"
                if os.path.exists(stats_path):
                    with open(stats_path) as handle:
                        line = handle.readline().rstrip("\n")
                    if line:
                        reason = line
                skipped_rows.append((sample_name, reason))

        os.makedirs(os.path.dirname(output.kept), exist_ok=True)
        with open(output.kept, "w") as handle:
            handle.write("\n".join(kept) + ("\n" if kept else ""))

        with open(output.skipped, "w") as handle:
            handle.write("sample\treason\n")
            for sample_name, reason in skipped_rows:
                handle.write(f"{sample_name}\t{reason}\n")

        with open(log[0], "w") as handle:
            handle.write(
                f"Kept {len(kept)}/{len(sample)} samples after trimming\n"
            )
            for sample_name, reason in skipped_rows:
                handle.write(f"SKIP {sample_name}: {reason}\n")

        if not kept:
            raise ValueError(
                "ERROR: No samples left after trimming. "
                f"See {output.skipped}"
            )


def get_kept_samples(wildcards=None):
    checkpoint_output = checkpoints.select_trimmed_samples.get()
    kept_path = checkpoint_output.output.kept
    with open(kept_path) as handle:
        return [line.strip() for line in handle if line.strip()]


rule fastqc_trimmed:
    input:
        fastq = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R{read}_trimmed.fastq.gz"),
        flag = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.done"),
        status = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.status"),
    output:
        html = os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.html"),
        zip = os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.zip")
    threads: 2
    resources:
        mem_mb=2000,
        runtime=60,
        nodes=1,
        cpus_per_task=2,
    log:
        os.path.join("logs", "fastqc_trimmed", "{sample}_R{read}.log")
    conda:
        "../envs/fastqc.yml"
    shell:
        """
        mkdir -p {work_dir}/metrics/trimmed

        if [ "$(head -n 1 {input.status})" != "ok" ]; then
            echo "SKIP: {wildcards.sample} was not trimmed successfully" | tee -a {log}
            exit 1
        fi

        if ! gunzip -t {input.fastq} 2>> {log}; then
            echo "ERROR: Input file {input.fastq} is corrupted" | tee -a {log}
            exit 1
        fi

        fastqc -t {threads} -o {work_dir}/metrics/trimmed {input.fastq} 2> {log}
        """
