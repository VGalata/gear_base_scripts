# About


# Structure

Custom processing/analysis/plotting scripts used for "XX".

## `R/`

Contains code of the R-package [superheat](https://github.com/rlbarter/superheat) which was slightly modified to accept additional arguments such as legend title and linkage method for `hclust`.

## `Scripts/`

- MIC value processing
    - `AST_process/10k_bacs_MIC_processing_v2.R`
- Parsing Kraken reports
    - `Taxonomy/parse_kraken_report.py`
- NCBI taxonomy from taxa assigned by Kraken
    - `Utils/get_full_tax.R`
- Filter samples by taxonomy and assembly quality
    - `Utils/filter_samples.R`
- Pan-genomes
    - Centroid HMMs
        - `Pan_genome/train_HMMs.py`
    - Centroid rates simulation
        - `Pan_genome/estimate_gene_rates.py`
    - Further analyses
        - Create `*.rds` data for scripts below
            - `SubProject_Pangenomes/scripts/create_meta.R`
            - `SubProject_Pangenomes/scripts/create_res_profiles.R`
        - Some utils
            - `SubProject_Pangenomes/scripts/utils.R`
        - Krona plot
            - `SubProject_Pangenomes/scripts/plot_tax_krona.R`
        - Processed resfams/essential genes mapping to centroids
            - `SubProject_Pangenomes/scripts/map_centroids2hits.py`
        - Processing essential gene hits
            - `SubProject_Pangenomes/scripts/process_essential_gene_hits.R`
        - Plotting essential gene hits
            - `SubProject_Pangenomes/scripts/plot_essential_gene_hits.R`
        - MSA of essential genes
            - `SubProject_Pangenomes/scripts/msa_essential_genes.py`
            - `SubProject_Pangenomes/scripts/cat_msa.py`
        - Process resfams hits
            - `SubProject_Pangenomes/scripts/process_resfams_hits.R`
        - Plot resfams hits
            - `SubProject_Pangenomes/scripts/plot_resfams_hits.R`
        - Plot MLST summary
            - `SubProject_Pangenomes/scripts/plot_mlst_sum.R`
        - Assembly quality plot
            - `SubProject_Pangenomes/scripts/plot_assembly_quali.R`
        - Centroid rates plot
            - `SubProject_Pangenomes/scripts/plot_pg_centroid_rates.R`
        - Pan-genome size overview
            - `SubProject_Pangenomes/scripts/plot_pg_size_ov.R`
        - Plot RaXML tree
            - `SubProject_Pangenomes/scripts/plot_raxml.R`
- Stat. analysis w/o genomic data
    - Res. vs. date (WMW tests)
        - `SubProject_ClassStat/Pipeline/scripts/data_analysis/res_date_tests.R`
    - Plot results
        - `SubProject_ClassStat/Pipeline/scripts/data_analysis/res_date_tests_plots.R`
    - Drug correlations (w/ plot)
        - `SubProject_ClassStat/Pipeline/scripts/data_analysis/drug_cor.R`
- Resistance associations
    - `Prediction/`
        - *Only EIGENSTRAT results on whole dataset were used*
