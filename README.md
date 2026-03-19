
# ADAMTS4-Stromal-Remodeling


Analyis code for: **Persistence of alveolar fibroblast-derived ADAMTS4+ cells in a preclinical model of delayed pulmonary fibrosis resolution** (Zabihi & Khadim et al., 2026)


### Project Overview

We identify a reversible switch between lipofibroblasts and myofibroblasts during pulmonary fibrosis. Using lineage tracing, scRNA-seq, and spatial transcriptomics, we show that ADAMTS4+ cells persist during delayed resolution and represent a therapeutic target in interstitial lung disease.


### Computational Workflow

We used a hybrid **Python/R** pipeline to analyze the single-cell and spatial data:


* **Python (Scanpy, scVelo, Squidpy)**: Preprocessing, sample integration, RNA velocity, and spatial mapping.

    * *References*: [Wolf et al. 2018](https://doi.org/10.1186/s13059-017-1382-0); [Bergen et al. 2020](https://doi.org/10.1038/s41587-020-0591-3); [Palla et al. 2022](https://doi.org/10.1038/s41592-021-01358-2)

* **R (Seurat, SeuratExtend, presto)**: Downstream visualization, stromal profiling, and fast differential expression.

    * *References*: [Butler et al. 2018](https://doi.org/10.1038/nbt.4096); [Hua et al. 2025](https://doi.org/10.1093/gigascience/giaf076); [Korsunsky et al. 2019](https://doi.org/10.1101/653253)

* **Statistics & Visualization**: GSEA via **clusterProfiler** and **fgsea**, with patterns visualized via **ComplexHeatmap** and **tidyplots**.

    * *References*: [Yu et al. 2012](https://doi.org/10.1089/omi.2011.0118); [Korotkevich et al. 2021](https://doi.org/10.1101/060012); [Gu et al. 2016](https://doi.org/10.1093/bioinformatics/btw313); [Engler 2025](https://doi.org/10.1002/imt2.70018)

---

### Analysis Pipeline

* **[notebooks/0-4]**: Murine scRNA-seq (Figs 2, S2–S4). Preprocessing, annotation, and trajectory mapping of *Fgf10*+ tdTom cells.

* **[notebooks/5-6]**: HLCA Comparative Analysis (Fig 7). Stromal sub-clustering and enrichment analysis across health & disease.

    * *Reference*: [Sikkema et al. *Nat Med* 2023](https://www.nature.com/articles/s41591-023-02327-2)

* **[notebooks/7]**: Spatial Transcriptomics (Fig 8). Visium analysis of human IPF tissues and *ADAMTS4/CTHRC1* mapping.

    * *Reference*: [Franzén et al. *Nat Genet* 2024](https://www.nature.com/articles/s41588-024-01819-2)

---

### Data

* **Data**: Both raw sequencing files and processed data objects are publicly available through the Gene Expression Omnibus (GEO) under accession number **[GSE295566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295566)**.


