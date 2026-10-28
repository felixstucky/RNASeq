rule fastqc:
    input:
        read1 = os.path.join(work_dir, symlink_dir, "{sample}_R1.fastq.gz"),
        read2 = os.path.join(work_dir, symlink_dir, "{sample}_R2.fastq.gz")
    params:
        outdir = os.path.join(work_dir, "metrics")
    output:
        html = [
            os.path.join(work_dir, "metrics", "{sample}_R1_fastqc.html"),
            os.path.join(work_dir, "metrics", "{sample}_R2_fastqc.html")
        ],
        zip = [
            os.path.join(work_dir, "metrics", "{sample}_R1_fastqc.zip"),
            os.path.join(work_dir, "metrics", "{sample}_R2_fastqc.zip")
        ]
    conda:
        "../envs/fastqc.yml"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=60,
        nodes=1,
        cpus_per_task=2,
    log:
        os.path.join(log_dir, "fastqc_{sample}.log")
    shell:
        """
        mkdir -p {params.outdir}
        echo "Starting FastQC for {wildcards.sample}..." | tee -a {log}
        fastqc --threads {threads} {input.read1} {input.read2} -o {params.outdir} 2>> {log}
        echo "✓ FastQC completed for {wildcards.sample}" | tee -a {log}
        """