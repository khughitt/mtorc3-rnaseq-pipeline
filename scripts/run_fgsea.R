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
  ungroup()

# convert to a named vector
gene_scores <- setNames(deseq2_res$log2FoldChange, deseq2_res$symbol)

# compute functional enrichment and store result
res <- fgsea(pathways = gene_sets, stats = gene_scores, eps = 0, nPermSimple = 10000) %>%
  arrange(padj)

write_tsv(res, snakemake@output[[1]])
