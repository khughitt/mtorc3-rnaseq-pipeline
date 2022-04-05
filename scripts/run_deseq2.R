#
# mTORC3 RNA-Seq Differential Expression Analysis (DESeq2)
#
library(annotables)
library(DESeq2)
library(tidyverse)

set.seed(1)

# load count table & sample metadata
sample_mdata <- read_tsv(snakemake@input[[2]], col_types='fff')

# add a combined ko/nutrient status variable
sample_mdata$condition <- factor(paste(sample_mdata$knockout, sample_mdata$nutrient, sep="_"))

raw_counts <- read_tsv(snakemake@input[[1]], col_types=cols()) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# remove genes with no reads mapped to them; this reduces memory required + they will be
# filtered out by deseq2 anyway
num_reads <- apply(raw_counts, 1, sum)

mask <- num_reads > 0
raw_counts <- raw_counts[mask, ]

# 1) effect of knockout/nutrient status, controlling for the other variable + batch
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = sample_mdata,
                              design= ~ batch + knockout + nutrient)
dds <- DESeq(dds)

# vectors to store contrast names and number of significant genes associated with each
# contrast; used to construct a summary table
contrasts1 <- c()
num_sig05 <- c()
num_sig01 <- c()
num_sig001 <- c()

# map ensgene ids to gene symbols (1:many)
gene_symbols <- grch38$symbol[match(rownames(raw_counts), grch38$ensgene)]

# iterate over constrasts and store each
out_dir <- dirname(snakemake@output[[1]])

ko_contrasts <- combn(unique(sample_mdata$knockout), 2)

for (i in 1:ncol(ko_contrasts)) {
  cond1 <- as.character(ko_contrasts[1, i])
  cond2 <- as.character(ko_contrasts[2, i])

  print(sprintf("Computing shrunken LFC for %s vs. %s...", cond1, cond2))

  # pos foldchange: lhs > rhs
  res <- lfcShrink(dds, contrast=c("knockout", cond1, cond2), type="ashr") %>%
    as.data.frame() %>%
    rownames_to_column("ensgene")

  res <- cbind(ensgene = res[, 1], symbol = gene_symbols, res[, 2:ncol(res)])

  num_sig05 <- c(num_sig05, sum(res$padj < 0.05, na.rm=TRUE))
  num_sig01 <- c(num_sig01, sum(res$padj < 0.01, na.rm=TRUE))
  num_sig001 <- c(num_sig001, sum(res$padj < 0.001, na.rm=TRUE))

  contrast <- sprintf("knockout_%s_vs_%s", cond2, cond1)
  contrasts1 <- c(contrasts1, contrast)

  write_tsv(res, file.path(out_dir, sprintf("%s.tsv", contrast)))
}

# starved vs. nutrient rich
res <- lfcShrink(dds, coef="nutrient_nutri_vs_starved", type="ashr") %>%
  as.data.frame() %>%
  rownames_to_column("ensgene")

res <- cbind(ensgene = res[, 1], symbol = gene_symbols, res[, 2:ncol(res)])

num_sig05 <- c(num_sig05, sum(res$padj < 0.05, na.rm=TRUE))
num_sig01 <- c(num_sig01, sum(res$padj < 0.01, na.rm=TRUE))
num_sig001 <- c(num_sig001, sum(res$padj < 0.001, na.rm=TRUE))

contrasts1 <- c(contrasts1, "nutrient_nutri_vs_starved")

write_tsv(res, file.path(out_dir, "nutrient_nutri_vs_starved.tsv"))

# 2) effect of knockout/nutrient status (single factor / combined model)
dds2 <- DESeqDataSetFromMatrix(countData = raw_counts,
                               colData = sample_mdata,
                               design= ~ batch + condition)
dds2 <- DESeq(dds2)

# contrasts of interest (combined conditions)
contrasts2 <- c("condition_m7_starved_vs_m7_nutri", "condition_Rap_nutri_vs_m7_nutri",
                "condition_Rap_starved_vs_m7_nutri", "condition_Ric_nutri_vs_m7_nutri",
                "condition_Ric_starved_vs_m7_nutri", "condition_v_nutri_vs_m7_nutri",
                "condition_v_starved_vs_m7_nutri")

# iterate over constrasts and store each
for (contrast in contrasts2) {
  print(sprintf("Computing shrunken LFC for %s...", contrast))

  res2 <- lfcShrink(dds2, coef=contrast, type="apeglm") %>%
    as.data.frame() %>%
    rownames_to_column("ensgene")

  res2 <- cbind(ensgene = res2[, 1], symbol = gene_symbols, res2[, 2:ncol(res2)])

  num_sig05 <- c(num_sig05, sum(res2$padj < 0.05, na.rm=TRUE))
  num_sig01 <- c(num_sig01, sum(res2$padj < 0.01, na.rm=TRUE))
  num_sig001 <- c(num_sig001, sum(res2$padj < 0.001, na.rm=TRUE))

  write_tsv(res2, file.path(out_dir, sprintf("%s.tsv", contrast)))
}

# generate & store summary table
summary <- data.frame(
  contrast = c(contrasts1, contrasts2),
  num_sig05,
  num_sig01,
  num_sig001
)

ind <- length(snakemake@output)
write_tsv(summary, snakemake@output[[ind]])
