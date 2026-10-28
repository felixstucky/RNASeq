rule create_links:
    input:
        metadata = metadata_file,
    output:
        r1 = expand(os.path.join(work_dir, symlink_dir, "{sample}_R1.fastq.gz"), sample=sample),
        r2 = expand(os.path.join(work_dir, symlink_dir, "{sample}_R2.fastq.gz"), sample=sample),
    log:
        os.path.join(work_dir, log_dir, "symlink.log")
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
        cpus_per_task=1,
    threads: 1
    run:
        import pandas as pd
        import os
        from pathlib import Path

        print(f"[LOG] Loading metadata from: {input.metadata}")
        df_metadata = pd.read_csv(input.metadata)

        required_cols = ["Sample_Name", "Read", "Absolute_Path"]
        missing_cols = [col for col in required_cols if col not in df_metadata.columns]
        if missing_cols:
            raise ValueError(f"ERROR: Metadata file missing required columns: {missing_cols}")

        sample_reads = df_metadata.groupby("Sample_Name")["Read"].apply(list)
        for sample_name, reads in sample_reads.items():
            if len(reads) != len(set(reads)):
                raise ValueError(f"ERROR: Duplicate reads for sample {sample_name}: {reads}")

        os.makedirs(os.path.join(work_dir, symlink_dir), exist_ok=True)

        successful_links = 0
        for index, row in df_metadata.iterrows():
            fastq_file = Path(row["Absolute_Path"])
            if not fastq_file.exists():
                raise FileNotFoundError(f"ERROR: Input file not found: {fastq_file}")

            output_file_Read = os.path.join(work_dir, symlink_dir, f'{row["Sample_Name"]}_{row["Read"]}.fastq.gz')

            if os.path.exists(output_file_Read):
                os.remove(output_file_Read)

            try:
                os.symlink(fastq_file, output_file_Read)
                successful_links += 1
                print(f"[LOG] Created symlink for {row['Sample_Name']}_{row['Read']}")
            except Exception as e:
                raise RuntimeError(f"ERROR: Failed to create symlink for {row['Sample_Name']}_{row['Read']}: {e}")

        print(f"[LOG] ✓ Successfully created {successful_links} symlinks")
        created_files = os.listdir(os.path.join(work_dir, symlink_dir))
        print(f"[LOG] Files in symlink directory: {created_files}")