checkpoint salmon_mapping_qc:
    input:
        meta_info = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "aux_info", "meta_info.json"),
            sample=get_kept_samples(wildcards),
        ),
        lib_format = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "lib_format_counts.json"),
            sample=get_kept_samples(wildcards),
        ),
        kept = os.path.join(work_dir, "metrics", "kept_samples_after_trim.txt"),
    output:
        tsv = os.path.join(work_dir, "metrics", "salmon_mapping_qc_mqc.tsv"),
        kept_after_mapping = os.path.join(
            work_dir, "metrics", "kept_samples_after_mapping.txt"
        ),
        skipped_after_mapping = os.path.join(
            work_dir, "metrics", "skipped_after_mapping.tsv"
        ),
        skipped_mqc = os.path.join(
            work_dir, "metrics", "salmon_mapping_excluded_mqc.tsv"
        ),
    params:
        min_mapping_rate = (
            "null" if min_mapping_rate in (None, False) else min_mapping_rate
        ),
        salmon_root = os.path.join(work_dir, salmon_dir),
        script = os.path.join(workflow.basedir, "scripts", "salmon_mapping_qc.py"),
    log:
        os.path.join(work_dir, log_dir, "salmon_mapping_qc.log"),
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
        echo "[LOG] Starting salmon_mapping_qc $(date -Is)" > {log}

        manifest="$(dirname {output.tsv})/salmon_mapping_qc_manifest.tsv"
        : > "$manifest"
        while IFS= read -r sample || [ -n "$sample" ]; do
            sample=$(echo "$sample" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            [ -z "$sample" ] && continue
            meta="{params.salmon_root}/$sample/aux_info/meta_info.json"
            lib="{params.salmon_root}/$sample/lib_format_counts.json"
            printf '%s\t%s\t%s\n' "$sample" "$meta" "$lib" >> "$manifest"
        done < {input.kept}

        echo "[LOG] Manifest lines: $(wc -l < "$manifest")" >> {log}
        python3 {params.script} \
            --manifest "$manifest" \
            --output {output.tsv} \
            --kept-output {output.kept_after_mapping} \
            --skipped-output {output.skipped_after_mapping} \
            --skipped-mqc-output {output.skipped_mqc} \
            --min-mapping-rate {params.min_mapping_rate} \
            >> {log} 2>&1

        echo "[LOG] Done salmon_mapping_qc $(date -Is)" >> {log}
        """


def get_mapping_passed_samples(wildcards=None):
    checkpoint_output = checkpoints.salmon_mapping_qc.get()
    kept_path = checkpoint_output.output.kept_after_mapping
    with open(kept_path) as handle:
        passed = [line.strip() for line in handle if line.strip()]
    if not passed:
        raise ValueError(
            "ERROR: No samples left after mapping QC. "
            f"See {checkpoint_output.output.skipped_after_mapping}"
        )
    return passed
