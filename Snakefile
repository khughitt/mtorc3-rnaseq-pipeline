"""
mTORC3 RNA-Seq Analysis
Mar 16, 2022

KH: Re-analyzing Joe's mTORC1/3 knockout data, starting from the HiSat2 mapped reads.

Outline:

- count tables (featureCounts)
- differential expression (DESeq2)
- gene set enrichment analysis (fgsea?)

TODO: create data packages for each output!
"""
import os
import pandas as pd
import yaml
from nodes import BioMatrixNode

# wildcards
knockouts = ['m7_d5', 'Rap_d1', 'Ric_d5', 'v_d4']
nutrients = ['starved', 'nutri']
batches = ['1', '2', '3', '4']

# msigdb
msigdb_gmt = "/data/raw/msigdb/v7.5.1/msigdb.v7.5.1.symbols.gmt"

# deseq contrasts
contrasts = ['condition_m7_starved_vs_m7_nutri',
             'condition_Rap_nutri_vs_m7_nutri', 'condition_Rap_starved_vs_m7_nutri',
             'condition_Ric_nutri_vs_m7_nutri', 'condition_Ric_starved_vs_m7_nutri',
             'condition_v_nutri_vs_m7_nutri', 'condition_v_starved_vs_m7_nutri',
             'knockout_Rap_vs_m7', 'knockout_Ric_vs_m7', 'knockout_v_vs_m7',
             'knockout_Ric_vs_Rap', 'knockout_v_vs_Rap', 'knockout_v_vs_Ric',
             'nutrient_nutri_vs_starved']

rule all:
    input:
        "/data/packages/mtorc3/raw-counts/datapackage.yml",
        "/data/packages/mtorc3/deseq2/summary.tsv",
        expand("/data/packages/mtorc3/fgsea/{contrast}.tsv", contrast=contrasts)

rule fgsea:
    input:
        "/data/packages/mtorc3/deseq2/{contrast}.tsv",
        "/data/packages/mtorc3/raw-counts/column-metadata.tsv"
    output:
        "/data/packages/mtorc3/fgsea/{contrast}.tsv"
    params:
        gmt=msigdb_gmt
    script:
        "scripts/run_fgsea.R"

rule deseq2:
    input:
        "/data/packages/mtorc3/raw-counts/data.tsv",
        "/data/packages/mtorc3/raw-counts/column-metadata.tsv",
    output:
        expand("/data/packages/mtorc3/deseq2/{contrast}.tsv", contrast=contrasts),
        "/data/packages/mtorc3/deseq2/summary.tsv"
    script:
        "scripts/run_deseq2.R"

rule package_raw_counts:
    input:
        "/data/packages/mtorc3/raw-counts/data.tsv",
        "/data/packages/mtorc3/raw-counts/column-metadata.tsv",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-library-sizes.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-pca.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-umap.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-heatmap-pearson.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-heatmap-spearman.png"
    output:
        "/data/packages/mtorc3/raw-counts/datapackage.yml"
    run:

        # load package metadata
        with open("metadata/metadata.yml") as fp:
            metadata = yaml.load(fp, Loader=yaml.FullLoader)

        # create node & export to data package
        node = BioMatrixNode(input[0], column_metadata=input[1], **metadata)
        node.to_pkg(os.path.dirname(output[0]), data_format="tsv")

rule raw_count_figures:
    input:
        "/data/packages/mtorc3/raw-counts/data.tsv",
        "/data/packages/mtorc3/raw-counts/column-metadata.tsv"
    output:
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-library-sizes.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-pca.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-umap.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-heatmap-pearson.png",
        "/data/packages/mtorc3/raw-counts/fig/raw-counts-sample-heatmap-spearman.png"
    script:
        "scripts/plot_raw_counts.R"

rule combine_counts:
    input:
        expand("/data/proj/mtorc3/subread/{knockout}_{nutrient}_{batch}_all_trimmed.txt",
               knockout=knockouts, nutrient=nutrients, batch=batches)
    output:
        "/data/packages/mtorc3/raw-counts/data.tsv",
        "/data/packages/mtorc3/raw-counts/column-metadata.tsv"
    run:
        dfs = [pd.read_csv(x, sep='\t', skiprows=1) for x in input]

        combined = pd.concat([dfs[0].Geneid] + [x.iloc[:, -1] for x in dfs], axis=1)
        combined.columns = [os.path.basename(x).replace("_all_trimmed.bam", "") for x in combined.columns]
        combined.to_csv(output[0], sep="\t", index=False)

        mdata_rows = []

        for sample_name in combined.columns[1:]:
            parts = sample_name.split('_')

            mdata_rows.append({
                "knockout": parts[0],
                "nutrient": parts[2],
                "batch": parts[3]
            })

        sample_metadata = pd.DataFrame.from_dict(mdata_rows)
        sample_metadata.to_csv(output[1], sep="\t", index=False)

rule count_reads:
    input:
        "/data/proj/mtorc3/hisat2/{knockout}_{nutrient}_{batch}_all_trimmed.bam"
    output:
        "/data/proj/mtorc3/subread/{knockout}_{nutrient}_{batch}_all_trimmed.txt"
    shell:
        """
        featureCounts \
            -p \
            -T {threads} \
            -t exon \
            -g gene_id \
            -a /data/ref/ensembl/GRCh38/105/Homo_sapiens.GRCh38.105.gtf \
            -o {output} \
            {input}
        """


# rule raw_counts_summary:
#     input:
#         "/data/packages/mtorc3/raw-counts/data.tsv",
#         "/data/packages/mtorc3/raw-counts/column-metadata.tsv"
#     output:
#         "/data/packages/mtorc3/summary/raw_counts_summary.html"
#     script:
#         "scripts/raw_counts_summary.Rmd"
