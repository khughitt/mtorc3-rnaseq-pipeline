#
# mTORC3 Raw Counts Summary Views
#
library(tidyverse)
library(uwot)
library(NMF)
library(viridis)

set.seed(1)

raw_counts <- read_tsv(snakemake@input[[1]], col_types = cols()) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# load count table & sample metadata
sample_mdata <- read_tsv(snakemake@input[[2]], col_types = 'fff')

lib_sizes <- apply(raw_counts, 2, sum)

dat <- cbind(sample_mdata, total_reads = lib_sizes)

ggplot(dat, aes(x = knockout, y = total_reads, fill = knockout, label = batch, group = batch)) +
  geom_col(position = position_dodge(0.9)) +
  facet_wrap(~nutrient) +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5) +
  theme_bw() +
  ggtitle("Library sizes")

ggsave(snakemake@output[[1]], width=1080, height=1080, units='px', dpi=192)

# remove genes with zero variance
num_reads <- apply(raw_counts, 1, sum)
raw_counts <- raw_counts[num_reads > 0, ]

pca <- prcomp(t(raw_counts), scale = TRUE)

pca_dat <- pca$x[, 1:2]
colnames(pca_dat) <- c("PC1", "PC2")

# compute variance explained
var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 2)

# add metadata
pca_dat <- pca_dat %>%
  as.data.frame() %>%
  rownames_to_column('row_label')

pca_dat <- cbind(pca_dat, sample_mdata)

ggplot(pca_dat, aes(x = PC1, y = PC2, color = knockout, shape = nutrient, label = batch)) +
  geom_point() +
  geom_text(hjust = 0, vjust = 0, nudge_y = 0.05, nudge_x = 0.05, size = 2) + 
  theme_bw() +
  ggtitle("Sample PCA") +
  xlab(sprintf("PC1 (%.2f%% variance)", var_explained[1])) +
  ylab(sprintf("PC2 (%.2f%% variance)", var_explained[2]))

ggsave(snakemake@output[[2]], width=1080, height=1080, units='px', dpi=192)

# UMAP
umap_dat <- umap(t(raw_counts), n_neighbors = 15, n_components = 2,
                 init = 'spectral', scale = TRUE, min_dist = 0.01)
colnames(umap_dat) <- c("UMAP1", "UMAP2")

# add metadata
umap_dat <- umap_dat %>%
  as.data.frame() %>%
  rownames_to_column('row_label')

umap_dat <- cbind(umap_dat, sample_mdata)

ggplot(umap_dat, aes(x = UMAP1, y = UMAP2, color = knockout, shape = nutrient, label = batch)) +
  geom_point() +
  geom_text(hjust = 0, vjust = 0, nudge_y = 0.05, nudge_x = 0.05, size = 2) + 
  theme_bw() +
  ggtitle("Sample UMAP") +
  xlab("UMAP1") +
  ylab("UMAP2")

ggsave(snakemake@output[[3]], width=1080, height=1080, units='px', dpi=192)

#
# Heatmap (pearson)
#
cor_mat <- cor(raw_counts)

aheatmap(cor_mat,
         annRow = sample_mdata[,1:2], annCol = sample_mdata[, 3],
         color = viridis(100), main = "Sample Correlation Heatmap (Pearson)",
         filename=snakemake@output[[4]])

# store sample correlation matrix
out_dir <- dirname(snakemake@output[[1]])

cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  write_tsv(file.path(out_dir, "pearson_cor_mat.tsv"))

cor_mat <- cor(raw_counts, method='spearman')

aheatmap(cor_mat,
         annRow = sample_mdata[,1:2], annCol = sample_mdata[, 3],
         color = viridis(100), main = "Sample Correlation Heatmap (Spearman)",
         filename=snakemake@output[[5]])

# store sample correlation matrix
cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  write_tsv(file.path(out_dir, "spearman_cor_mat.tsv"))
