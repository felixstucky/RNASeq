args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  stop("Expected 8 arguments: quant_files, sample_names, tx2gene, transcript_counts_out, transcript_tpm_out, gene_counts_out, gene_counts_non_zero_out, gene_counts_pc_out")
}

quant_files <- strsplit(args[[1]], ",", fixed = TRUE)[[1]]
sample_names <- strsplit(args[[2]], ",", fixed = TRUE)[[1]]
tx2gene_file <- args[[3]]
transcript_counts_out <- args[[4]]
transcript_tpm_out <- args[[5]]
gene_counts_out <- args[[6]]
gene_counts_non_zero_out <- args[[7]]
gene_counts_pc_out <- args[[8]]

suppressPackageStartupMessages(library(tximport))

if (length(quant_files) != length(sample_names)) {
  stop("quant_files and sample_names must have the same length")
}

names(quant_files) <- sample_names
tx2gene <- read.delim(tx2gene_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("transcript_id", "gene_id", "gene_name", "gene_type")
missing_cols <- setdiff(required_cols, colnames(tx2gene))
if (length(missing_cols) > 0) {
  stop(paste("Missing required tx2gene columns:", paste(missing_cols, collapse = ", ")))
}

tx2gene <- tx2gene[!duplicated(tx2gene$transcript_id), ]

txi_tx <- tximport(quant_files, type = "salmon", txOut = TRUE, countsFromAbundance = "no")
transcript_counts <- data.frame(transcript_id = rownames(txi_tx$counts), txi_tx$counts, check.names = FALSE)
transcript_tpm <- data.frame(transcript_id = rownames(txi_tx$abundance), txi_tx$abundance, check.names = FALSE)

tx2gene_min <- tx2gene[, c("transcript_id", "gene_id")]
txi_gene <- tximport(quant_files, type = "salmon", tx2gene = tx2gene_min, countsFromAbundance = "no")

gene_counts <- data.frame(gene_id = rownames(txi_gene$counts), txi_gene$counts, check.names = FALSE)
gene_meta <- unique(tx2gene[, c("gene_id", "gene_name", "gene_type")])
gene_counts <- merge(gene_meta, gene_counts, by = "gene_id", all.y = TRUE, sort = FALSE)
gene_counts <- gene_counts[, c("gene_id", "gene_name", "gene_type", sample_names)]

count_only <- gene_counts[, sample_names, drop = FALSE]
non_zero <- rowSums(count_only) > 0
gene_counts_non_zero <- gene_counts[non_zero, , drop = FALSE]
gene_counts_pc <- gene_counts_non_zero[gene_counts_non_zero$gene_type == "protein_coding", , drop = FALSE]

dir.create(dirname(transcript_counts_out), recursive = TRUE, showWarnings = FALSE)
write.csv(transcript_counts, transcript_counts_out, row.names = FALSE, quote = FALSE)
write.csv(transcript_tpm, transcript_tpm_out, row.names = FALSE, quote = FALSE)
write.csv(gene_counts, gene_counts_out, row.names = FALSE, quote = FALSE)
write.csv(gene_counts_non_zero, gene_counts_non_zero_out, row.names = FALSE, quote = FALSE)
write.csv(gene_counts_pc, gene_counts_pc_out, row.names = FALSE, quote = FALSE)