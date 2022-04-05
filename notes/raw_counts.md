Raw Counts
==========

Possible Raptor Outliers
------------------------

Based on the sample correlation heatmap and PCA/UMAP plots for the raw counts, there appear
to be 1-2 _Raptor_ knockout samples which behave somewhat differently than the rest:

1. Rap_d1_starved_2
2. Rap_d1_nutri_3

Overall, the mean sample correlations are still fairly similar.

Average Pearson correlation for each Raptor sample with all other Raptor samples:

```
> colMeans(rap_cor_mat)
Rap_d1_starved_1 Rap_d1_starved_2 Rap_d1_starved_3 Rap_d1_starved_4   Rap_d1_nutri_1   Rap_d1_nutri_2   Rap_d1_nutri_3   Rap_d1_nutri_4
       0.9827347        0.9504418        0.9850015        0.9839592        0.9808928        0.9781916        0.9327462        0.9810346
```

Interestingly, only the `Rap_d1_starved_2` sample really stands out in the PCA/UMAP plots.

The `Rap_d1_nutri_3` sample _is_ on the extreme end of PC2, in the PCA plot, although it
still groups with the other Raptor samples, overall.

The `Rap_d1_starved_2` does also have a somewhat lower total number of reads, compared
with the other Raptor samples.

Finally, when assessing the maximum number of reads mapped to a single gene, for each
sample, `Rap_d1_nutri_3` has a much higher value (~670,000) compared to the other
Raptor/nutrient samples (~410-460,000).
