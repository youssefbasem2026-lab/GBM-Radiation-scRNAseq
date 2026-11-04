<img width="1310" height="532" alt="image" src="https://github.com/user-attachments/assets/91984392-5385-486c-9e08-0eb347b6c41e" /># Single-Cell Transcriptomic Reanalysis of Radiation Response in Glioblastoma

**Dataset:** GEO [GSE162931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162931)  
**Analysis:** GBM131 and GBM022, untreated (Naive) versus two days after radiation (RT2d)  
**Author:** [Your full name]  
**GitHub:** [youssefbasem2026-lab](https://github.com/youssefbasem2026-lab)

## Project Description

This repository documents an exploratory reanalysis of publicly available single-cell RNA sequencing (scRNA-seq) data from cultured glioblastoma (GBM) cells before and after radiation treatment. The workflow examines cellular transcriptional heterogeneity, GBM-associated expression programs, radiation-response genes, changes in cluster composition, and predicted ligand–receptor signaling using R-based single-cell tools.

The project is a **secondary analysis of public data**. It does not claim to have generated the original sequencing dataset or experimentally validated the computational predictions.

## Background

Glioblastoma is a highly heterogeneous malignant brain tumor. Radiation is a major component of treatment, but transcriptional responses may differ between malignant cell states and patients. Single-cell RNA sequencing can describe this heterogeneity and generate hypotheses about radiation-associated changes that may be obscured in bulk measurements.

## Rationale and Research Questions

The analysis asks:

1. How are transcriptional clusters distributed across the two patient-derived GBM cultures?
2. How do observed cluster proportions and GBM-associated expression programs differ between untreated and irradiated samples?
3. Which genes and pathways show candidate radiation-associated expression changes in a paired, sample-level analysis?
4. Which ligand–receptor interactions are computationally predicted among transcriptionally defined clusters?

These are **exploratory questions**; the study is not powered for broad clinical inference.

## Data Availability and Sample Description

The source data are publicly available through the NCBI Gene Expression Omnibus under accession **GSE162931**.

| Sample | GEO accession | Patient-derived culture | Condition |
|---|---|---|---|
| GBM131_Naive | GSM4967236 | GBM131 | Untreated |
| GBM131_RT2d | GSM4967237 | GBM131 | 2 days after radiation |
| GBM022_Naive | GSM4967242 | GBM022 | Untreated |
| GBM022_RT2d | GSM4967244 | GBM022 | 2 days after radiation |

Only these four samples are included in this workflow. Please consult the GEO record and associated original publication for experimental details and the correct citation of the dataset creators.

## Methods

### 1. Data Import and Quality Control

Count matrices are loaded into Seurat independently for each sample. The script calculates mitochondrial transcript percentages (`^MT-`) and examines gene/UMI distributions. Cells are filtered using **nFeature_RNA > 200**, **nFeature_RNA < 6000**, and **percent.mt < 10**. These are fixed exploratory thresholds rather than sample-adaptive thresholds.

### 2. Doublet Detection

Individual Seurat objects are converted to `SingleCellExperiment`, and `scDblFinder` is applied separately to each sample. Cells classified as singlets are retained. Pre- and post-filtering plots are generated for inspection.

### 3. Normalization, Feature Selection, and Dimensional Reduction

The workflow uses Seurat normalization, selection of 2,000 highly variable features, scaling, and principal component analysis. Sample-level exploratory clustering and UMAP visualization precede the merged analysis.

### 4. Merging and Harmony Integration

The four samples are merged with patient and treatment metadata. Harmony is applied using **patient identity** as the grouping variable, followed by graph-based clustering and UMAP embedding using the first 20 Harmony dimensions. Because patient identity and tumor biology may be linked, correction should be inspected for possible over-integration.

### 5. Cluster Annotation and GBM State Markers

SingleR with the Human Primary Cell Atlas reference provides provisional reference labels. In cultured malignant GBM cells, these labels **must not be interpreted as validated normal cell identities**. Marker expression and exploratory OPC-like, NPC-like, AC-like, and MES-like gene programs are examined using dot plots, feature plots, violin plots, and module scores. The short marker lists are **proxies, not validated full GBM state signatures**.

### 6. Differential Abundance and Cluster Composition

MiloR is used for neighborhood-level exploratory differential abundance, and stacked bar plots summarize cluster composition across samples. With only two paired patient-derived cultures, neighborhood significance estimates should be interpreted cautiously.

### 7. Paired Pseudobulk Differential Expression

Counts are aggregated by biological sample. edgeR models patient identity and treatment condition (`~ patient_id + condition`) to assess candidate RT-versus-Naive differences. Outputs include differential-expression tables, volcano plots, MA plots, MDS plots, and a pseudobulk expression heatmap.

**Statistical caution:** There are only **two paired biological patient samples**. Thousands of cells do not constitute thousands of independent biological replicates; statistical estimates and adjusted p-values may be unstable.

### 8. Functional Enrichment and Gene Set Analysis

Gene Ontology (Biological Process, Molecular Function, and Cellular Component), KEGG over-representation analysis, and ranked GO gene set enrichment are implemented with `clusterProfiler`, `org.Hs.eg.db`, and `enrichplot`. Enrichment is conditional on the differential-expression output and should be treated as hypothesis-generating.

### 9. Predicted Cell–Cell Communication

CellChat estimates potential ligand–receptor interactions between computationally defined clusters using the human interaction database. Network diagrams and interaction tables summarize predicted communication. Because these are **cultured GBM samples**, the results should not be described as an in vivo tumor microenvironment or proof of physical signaling between cell populations.

### 10. Radiation Response, Cell Cycle, and Transcription Factors

The script visualizes selected DNA damage/radiation-response genes, scores cell-cycle phases, and displays RNA expression of selected transcription factors. **TF expression is not equivalent to inferred TF activity.** Optional Slingshot pseudotime is exploratory and cannot demonstrate that radiation caused cell-state transitions.

## Results and Figures

> **To be populated from verified analysis outputs.** The analysis script defines the outputs below, but their successful generation and biological interpretation must be confirmed before reporting findings.

| Analysis output | Example output location | Status |
|---|---|---|
| Integrated UMAP and cluster visualization | `GBM_plots/` | Add verified figure |
| Cluster composition by sample | `GBM_plots/22_Cluster_composition.png` | Add verified figure |
| MiloR differential abundance | `GBM_plots/23_Milo_DA_scatter.png` | Add verified figure |
| Pseudobulk library sizes | `GBM_plots/24_Pseudobulk_library_sizes.png` | Add verified figure |
| RT vs Naive volcano | `GBM_plots/26_Volcano.png` | Add verified figure |
| Pseudobulk DE table | `GBM_RT_vs_Naive_DEG.csv` | Verify results |
| CellChat communication network | `GBM_plots/28_CellChat_Communication_Network.png` | Add verified figure |
| GBM state scores | `GBM_OPC_like_sample_scores.csv` and corresponding state files | Verify results |
| Radiation-response genes | `GBM_plots/36_Radiation_DotPlot.png` | Add verified figure |
| Cell-cycle composition | `GBM_plots/37_CellCycle_Composition.png` | Add verified figure |
| GO gene set enrichment | `GBM_plots/38_GSEA_GO_DotPlot.png` | Add verified figure |

### Suggested Figure Gallery

Once the actual PNG files are uploaded, uncomment the following examples:

<!--
![Cluster composition by sample](GBM_plots/22_Cluster_composition.png)
![RT versus Naive volcano plot](GBM_plots/26_Volcano.png)
![Radiation response markers](GBM_plots/36_Radiation_DotPlot.png)
-->

## Discussion

The analysis is designed to investigate heterogeneity in cultured GBM cells and candidate transcriptional changes following irradiation. Any observed differences must be evaluated against patient-specific variation, cluster composition, sequencing quality, and the very small number of independent samples. Pathway enrichment and ligand–receptor predictions identify hypotheses rather than mechanisms established by experiment.

**Interpretation of specific genes, pathways, and treatment effects will be added after reviewing the actual result tables and plots.**

## Limitations

- **Two patient pairs:** The sample size limits statistical power, uncertainty estimation, and generalizability.
- **Cultured cells:** The dataset does not reconstruct the immune or stromal architecture of intact human tumors.
- **Reference annotation:** SingleR labels from normal-cell references are provisional for malignant GBM cells.
- **Integration effects:** Harmony correction may remove patient-associated biological differences and requires diagnostic inspection.
- **CellChat predictions:** Transcript-based ligand–receptor scores do not establish functional or physical communication.
- **Gene program scoring:** Small curated marker sets are not substitutes for validated comprehensive GBM-state signatures.
- **No spatial coordinates:** UMAP is a transcriptional embedding, **not spatial transcriptomics**.
- **Pseudotime:** Computational ordering does not establish temporal or causal transitions after radiation.
- **Reproducibility:** Package versions, input file layout, and all script sections should be checked on a clean R session before claiming a fully reproducible pipeline.

## Conclusion

This repository provides an R-based exploratory workflow for examining single-cell transcriptional heterogeneity and radiation-associated candidate responses in two patient-derived glioblastoma cultures. The outputs can guide further validation, but the current dataset is insufficient for definitive claims about radiation resistance mechanisms or clinical biomarkers.

## Repository Structure

```text
GBM-Radiation-scRNAseq/
├── README.md
├── GBM_scRNAseq_Analysis.R
├── requirements.R
├── .gitignore
└── GBM_plots/                # optional: selected verified figures
```

## How to Run

1. Install R and the packages used in the script (see `requirements.R`).
2. Download the relevant processed count matrices from GEO GSE162931.
3. Organize the four sample folders as expected by the script:

```text
data/extracted/
├── GBM131_Naive/
├── GBM131_RT2d/
├── GBM022_Naive/
└── GBM022_RT2d/
```

Each sample folder should contain the 10x-format matrix, barcodes, and feature files expected by `Read10X`.

4. Set the R working directory to the repository root, check file paths, and run the script in order. Review optional package requirements for CellChat and Slingshot.
5. Confirm each exported table and figure before using it in a report or publication.

> **Note:** This repository is an analysis record and may require code cleanup and version-specific adjustments before it can run end-to-end without intervention.

## Software and Packages

Core tools include **R, Seurat, SingleCellExperiment, scDblFinder, SingleR, celldex, Harmony, miloR, edgeR, clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2, patchwork, and CellChat**. Optional trajectory analysis uses **Slingshot**.

## Data Citation and Acknowledgments

The underlying sequencing data were generated by the investigators who deposited **GSE162931** in GEO. Please cite the dataset's associated primary publication and the relevant software papers when reusing the analysis. This repository represents an independent secondary computational reanalysis.

## Author

**Youssef Basem**  
GitHub: [@youssefbasem2026-lab](https://github.com/youssefbasem2026-lab)

## License

No license has been assigned yet. A code license can be selected after confirming the intended reuse permissions. The original dataset retains its own attribution and applicable terms.
