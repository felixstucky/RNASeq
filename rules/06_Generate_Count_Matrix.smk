rule generate_count_matrix:
    input:
        quant = lambda wildcards: expand(
            os.path.join(work_dir, salmon_dir, "{sample}", "quant.sf"),
            sample=get_mapping_passed_samples(wildcards),
        ),
        tx2gene = os.path.join(salmon_ref_dir, "tx2gene.tsv"),
        mapping_kept = os.path.join(work_dir, "metrics", "kept_samples_after_mapping.txt"),
        mapping_qc = os.path.join(work_dir, "metrics", "salmon_mapping_qc_mqc.tsv"),
    output:
        transcript_counts = os.path.join(work_dir, count_matrix_dir, "TranscriptCountMatrix.csv"),
        transcript_tpm = os.path.join(work_dir, count_matrix_dir, "TranscriptTPMMatrix.csv"),
        count_Matrix = os.path.join(work_dir, count_matrix_dir, "GeneCountMatrix.csv"),
        count_Matrix_non_zero = os.path.join(work_dir, count_matrix_dir, "GeneCountMatrix_non_zero_count.csv"),
        count_Matrix_protein_coding_non_zero = os.path.join(work_dir, count_matrix_dir, "GeneCountMatrix_non_zero_count_protein_coding_only.csv")
    params:
        work_dir = work_dir,
        count_matrix_dir = count_matrix_dir,
        quant_files = lambda wildcards, input: ",".join(input.quant),
        sample_names = lambda wildcards: ",".join(get_mapping_passed_samples(wildcards)),
    log:
        os.path.join(work_dir, log_dir, "gen_count_matrix.log")
    resources:
        mem_mb=16000,
        runtime=60,
        nodes=1,
        cpus_per_task=1,
    threads: 1
    conda:
        "../envs/r_tximport.yml"
    shell:
        """
        Rscript --vanilla {workflow.basedir}/scripts/salmon_tximport.R \
            "{params.quant_files}" \
            "{params.sample_names}" \
            {input.tx2gene} \
            {output.transcript_counts} \
            {output.transcript_tpm} \
            {output.count_Matrix} \
            {output.count_Matrix_non_zero} \
            {output.count_Matrix_protein_coding_non_zero} \
            > {log} 2>&1
        """
