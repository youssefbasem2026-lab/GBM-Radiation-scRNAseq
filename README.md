# Glioblastoma Single-Cell RNA-seq Analysis: Radiation Response

**Dataset:** [GSE162931 (NCBI GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162931)  
**Analysis author / repository maintainer:** **[YOUR NAME]**  
**Data:** Publicly available, independently generated dataset; credit the original dataset authors when using or publishing results.

## Overview

Exploratory R/Seurat workflow analyzing cultured glioblastoma (GBM) single-cell RNA-seq data before irradiation (Naive) and two days after irradiation (RT2d). The script includes four samples from two patient-derived models: GBM131 and GBM022. This repository contains an analysis workflow, **not a new spatial transcriptomics dataset**.

## Included samples

| Sample | Condition |
| --- | --- |
| GBM131_Naive | Untreated |
| GBM131_RT2d | Two days post-radiation |
| GBM022_Naive | Untreated |
| GBM022_RT2d | Two days post-radiation |

## Analysis workflow

1. Read 10X matrices, QC, cell filtering and scDblFinder doublet detection.
2. Normalize, select variable genes, PCA, clustering and UMAP.
3. Merge samples, integrate with Harmony and perform exploratory SingleR annotation.
4. Identify cluster markers and investigate GBM-associated gene programs.
5. Test neighborhood abundance with MiloR; perform paired pseudobulk differential expression with edgeR.
6. Visualize DE genes and perform GO/KEGG enrichment and GSEA.
7. Explore CellChat-inferred signaling, radiation-response genes, cell-cycle phases and optional pseudotime.

## Setup

Install the R and Bioconductor packages listed in [`requirements.R`](requirements.R). Package compatibility depends on your R/Bioconductor versions.

Download the relevant sample files from GEO and arrange extracted 10X files as:

```text
data/extracted/
  GBM131_Naive/
  GBM131_RT2d/
  GBM022_Naive/
  GBM022_RT2d/
```

Each folder should contain `matrix.mtx`, `barcodes.tsv`, and `features.tsv` or `genes.tsv` (possibly compressed). Run `GBM_scRNAseq_Analysis.R` **section by section** in RStudio, inspecting output and errors. The script assumes an output directory `GBM_plots/` exists before the first plot export; create it at the start.

```r
dir.create('GBM_plots', showWarnings = FALSE)
```

**Reproducibility status:** This is a cleaned publication copy of the user's exploratory script, **not an independently tested end-to-end pipeline**. Some sections may require debugging, package-specific adjustments, and manual review before reuse. Do not claim that all analyses completed successfully solely because they are included in the script.

## Important interpretation limits

- There are **two biological patient pairs**; differential-expression and abundance inferences are exploratory and statistical power is very limited.
- The source consists of **cultured GBM cells**, not spatially resolved tumor sections. UMAP is not a spatial map.
- SingleR labels based on normal reference atlases are not confirmed malignant GBM cell identities.
- CellChat estimates potential signaling from RNA expression; it does not establish experimentally observed interactions or immune-cell communication in this culture.
- Small GBM marker sets approximate cell-state programs and are not full validated signatures.
- Pseudotime does not demonstrate radiation-induced lineage transitions.

## Outputs

The workflow writes figures into `GBM_plots/` and various CSV and RDS analysis outputs. Large intermediate data, original raw counts, and RDS objects are intentionally excluded from version control.

## Attribution and citation

- **Analysis code and organization:** [YOUR NAME].
- **Original sequencing data:** Cite the original study associated with [GSE162931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162931). This repository does not claim ownership of those data.
- Cite the R packages and methods used in any manuscript.

## License

No license is assigned yet. Add a license only after confirming you own or have permission to distribute the analysis code, including any portions adapted from others.
