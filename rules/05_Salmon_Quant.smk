rule prepare_salmon_reference:
    input:
        genome = genome,
        gtf = gtf
    output:
        transcriptome = os.path.join(salmon_ref_dir, "transcriptome.fa"),
        decoys = os.path.join(salmon_ref_dir, "decoys.txt"),
        gentrome = os.path.join(salmon_ref_dir, "gentrome.fa"),
        tx2gene = os.path.join(salmon_ref_dir, "tx2gene.tsv")
    log:
        os.path.join(work_dir, log_dir, "prepare_salmon_reference.log")
    conda:
        "../envs/salmon.yml"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=240,
        nodes=1,
        cpus_per_task=2,
    shell:
        r"""
        mkdir -p {salmon_ref_dir}

        echo "[LOG] Building Salmon transcriptome reference..." | tee -a {log}

        gffread -w {output.transcriptome} -g {input.genome} {input.gtf} 2>> {log}

        grep '^>' {input.genome} | sed 's/^>//; s/ .*//' > {output.decoys}

        cat {output.transcriptome} {input.genome} > {output.gentrome}

        python - <<'PY'
from pathlib import Path

gtf_path = Path(r"{input.gtf}")
tx2gene_path = Path(r"{output.tx2gene}")

def parse_attributes(field):
    attributes = dict()
    for item in field.strip().strip(';').split(';'):
        item = item.strip()
        if not item or ' ' not in item:
            continue
        key, value = item.split(' ', 1)
        attributes[key] = value.strip().strip('"')
    return attributes

records = dict()
with gtf_path.open() as handle:
    for line in handle:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        feature = parts[2]
        if feature not in ('transcript', 'exon'):
            continue
        attributes = parse_attributes(parts[8])
        transcript_id = attributes.get('transcript_id')
        gene_id = attributes.get('gene_id')
        if not transcript_id or not gene_id:
            continue
        gene_name = attributes.get('gene_name', gene_id)
        gene_type = attributes.get('gene_biotype', attributes.get('gene_type', 'unknown'))
        records[transcript_id] = (gene_id, gene_name, gene_type)

with tx2gene_path.open('w') as handle:
    handle.write('transcript_id\tgene_id\tgene_name\tgene_type\n')
    for transcript_id in sorted(records):
        gene_id, gene_name, gene_type = records[transcript_id]
        handle.write(transcript_id + '\t' + gene_id + '\t' + gene_name + '\t' + gene_type + '\n')
PY

        echo "[LOG] Salmon reference prepared at {salmon_ref_dir}" | tee -a {log}
        """

rule salmon_index:
    input:
        gentrome = os.path.join(salmon_ref_dir, "gentrome.fa"),
        decoys = os.path.join(salmon_ref_dir, "decoys.txt")
    output:
        index = directory(salmon_index_dir)
    log:
        os.path.join(work_dir, log_dir, "salmon_index.log")
    conda:
        "../envs/salmon.yml"
    threads: 10
    resources:
        mem_mb=32000,
        runtime=240,
        nodes=1,
        cpus_per_task=10,
    shell:
        """
        mkdir -p {output.index}

        salmon index \
            -t {input.gentrome} \
            -d {input.decoys} \
            -i {output.index} \
            -k {salmon_kmer} \
            -p {threads} \
            &> {log}
        """

rule salmon_quant:
    input:
        index = salmon_index_dir,
        read1 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R1_trimmed.fastq.gz"),
        read2 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R2_trimmed.fastq.gz"),
        flag = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.done"),
        status = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.status"),
    output:
        quant = os.path.join(work_dir, salmon_dir, "{sample}", "quant.sf"),
        meta_info = os.path.join(work_dir, salmon_dir, "{sample}", "aux_info", "meta_info.json"),
        lib_format = os.path.join(work_dir, salmon_dir, "{sample}", "lib_format_counts.json")
    log:
        os.path.join(work_dir, log_dir, "salmon_{sample}.log")
    conda:
        "../envs/salmon.yml"
    threads: 8
    resources:
        mem_mb=64000,
        runtime=120,
        nodes=1,
        cpus_per_task=8,
    shell:
        """
        set -uo pipefail
        if [ "$(head -n 1 {input.status})" != "ok" ]; then
            echo "ERROR: refusing Salmon for skipped sample {wildcards.sample}" | tee -a {log}
            exit 1
        fi
        outdir=$(dirname "{output.quant}")
        mkdir -p "$outdir"
        echo "[LOG] Starting Salmon quant for {wildcards.sample}" | tee -a {log}
        
        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.read1} \
            -2 {input.read2} \
            -p {threads} \
            --validateMappings \
            --gcBias \
            --seqBias \
            --minAssignedFrags 1 \
            -o "$outdir" \
            >> {log} 2>&1
        salmon_status=$?
        
        if [ "$salmon_status" -ne 0 ]; then
            echo "ERROR: salmon quant failed for {wildcards.sample} with exit code $salmon_status" | tee -a {log}
            ls -lh "$outdir" >> {log} 2>&1 || true
            exit "$salmon_status"
        fi
        
        if [ ! -s {output.quant} ]; then
            echo "ERROR: Salmon quant.sf is empty or missing" | tee -a {log}
            ls -lh "$outdir" >> {log} 2>&1 || true
            exit 1
        fi

        if [ ! -s {output.meta_info} ]; then
            echo "ERROR: Salmon meta_info.json is empty or missing" | tee -a {log}
            ls -lh "$outdir" "$outdir/aux_info" >> {log} 2>&1 || true
            exit 1
        fi

        if [ ! -s {output.lib_format} ]; then
            echo "ERROR: Salmon lib_format_counts.json is empty or missing" | tee -a {log}
            ls -lh "$outdir" >> {log} 2>&1 || true
            exit 1
        fi
        """