#
# fgsea gene set enrichment
#
library(fgsea)
library(tidyverse)
library(GSEABase)

set.seed(321)

# load MSigDB gene sets
gene_sets <- geneIds(getGmt(snakemake@params$gmt))

# load deseq2 results
deseq2_res <- read_tsv(snakemake@input[[1]], col_types = cols())

# collapse multiple entries mapping to the same gene symbols;
# mostly impacts non-coding RNA genes (SNORDs, Y_RNA, etc.)
deseq2_res <- deseq2_res %>%
  group_by(symbol) %>%
  summarize(log2FoldChange=mean(log2FoldChange)) %>%
  arrange(log2FoldChange) %>%
  ungroup()

# drop entry with missing gene symbol
deseq2_res <- deseq2_res %>%
  filter(!is.na(symbol))

# convert to a named vector
gene_scores <- setNames(deseq2_res$log2FoldChange, deseq2_res$symbol)

# compute functional enrichment and store result
# res <- fgsea(pathways = gene_sets, stats = gene_scores, eps = 0, nPermSimple = 10000) %>%
#   arrange(padj)
res <- fgsea(pathways = gene_sets, stats = gene_scores, eps = 0, nPermSimple = 10000)

# in addition to the full results table, also create one with only significant &
# non-redundant gene sets included
res_sig <- res[order(pval)][padj < 0.01]

# collapsed <- collapsePathways(res_sig, gene_sets, gene_scores, nperm = 10000)
collapsed <- collapsePathways(res_sig, gene_sets, gene_scores)

write_tsv(res, snakemake@output[[1]])
