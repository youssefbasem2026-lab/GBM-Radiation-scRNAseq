# Single-Cell Transcriptomic Reanalysis of Radiation Response in Glioblastoma

**Dataset:** GEO [GSE162931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162931)  
**Analysis:** GBM131 and GBM022, untreated (Naive) versus two days after radiation (RT2d)  
**Author:** GitHub user [@youssefbasem2026-lab](https://github.com/youssefbasem2026-lab) (full name not provided)  
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

The following **63 figures** are the actual screenshots supplied for this reanalysis, preserved without altering plot contents. File names have been standardized for GitHub. Each panel is presented as an analysis output, not as proof of independent experimental validation. The order follows the workflow; where exact plot-to-sample correspondence is not visible, captions avoid assuming it.

### Quality Control — before and after filtering

The violin plots show the distributions of detected genes (nFeature_RNA), total UMI counts (nCount_RNA), and mitochondrial read fraction (percent.mt). Early panels include high-mitochondrial outliers; later panels show distributions after QC filtering. The sequence includes multiple intermediate views for GBM131 Naive.

**Figure 1. QC PreFilter GBM131 Naive.**

![Figure 1: QC PreFilter GBM131 Naive](figures/01_QC_PreFilter_GBM131_Naive.png)

**Figure 2. QC PreFilter GBM131 RT2d.**

![Figure 2: QC PreFilter GBM131 RT2d](figures/02_QC_PreFilter_GBM131_RT2d.png)

**Figure 3. QC PreFilter GBM022 Naive.**

![Figure 3: QC PreFilter GBM022 Naive](figures/03_QC_PreFilter_GBM022_Naive.png)

**Figure 4. QC PreFilter GBM022 RT2d.**

![Figure 4: QC PreFilter GBM022 RT2d](figures/04_QC_PreFilter_GBM022_RT2d.png)

**Figure 5. QC PostFilter GBM131 Naive Overview.**

![Figure 5: QC PostFilter GBM131 Naive Overview](figures/05_QC_PostFilter_GBM131_Naive_Overview.png)

**Figure 6. QC PostFilter GBM131 Naive Refined.**

![Figure 6: QC PostFilter GBM131 Naive Refined](figures/06_QC_PostFilter_GBM131_Naive_Refined.png)

**Figure 7. QC PostFilter GBM131 RT2d.**

![Figure 7: QC PostFilter GBM131 RT2d](figures/07_QC_PostFilter_GBM131_RT2d.png)

**Figure 8. QC PostFilter GBM022 Naive.**

![Figure 8: QC PostFilter GBM022 Naive](figures/08_QC_PostFilter_GBM022_Naive.png)

**Figure 9. QC PostFilter GBM022 RT2d.**

![Figure 9: QC PostFilter GBM022 RT2d](figures/09_QC_PostFilter_GBM022_RT2d.png)

### Doublet detection and removal

scDblFinder doublet scores and class assignments are visualized before and after removal. The retained-cell bar chart compares the four samples; bar heights are descriptive and are not a formal treatment-effect test.

**Figure 10. DoubletDetection GBM131 Naive Scatter.**

![Figure 10: DoubletDetection GBM131 Naive Scatter](figures/10_DoubletDetection_GBM131_Naive_Scatter.png)

**Figure 11. DoubletDetection GBM131 Naive Violin.**

![Figure 11: DoubletDetection GBM131 Naive Violin](figures/11_DoubletDetection_GBM131_Naive_Violin.png)

**Figure 12. DoubletDetection GBM131 RT2d Scatter.**

![Figure 12: DoubletDetection GBM131 RT2d Scatter](figures/12_DoubletDetection_GBM131_RT2d_Scatter.png)

**Figure 13. DoubletDetection GBM131 RT2d Violin.**

![Figure 13: DoubletDetection GBM131 RT2d Violin](figures/13_DoubletDetection_GBM131_RT2d_Violin.png)

**Figure 14. DoubletDetection GBM022 Naive Scatter.**

![Figure 14: DoubletDetection GBM022 Naive Scatter](figures/14_DoubletDetection_GBM022_Naive_Scatter.png)

**Figure 15. DoubletRemoval GBM131 Naive Scatter.**

![Figure 15: DoubletRemoval GBM131 Naive Scatter](figures/15_DoubletRemoval_GBM131_Naive_Scatter.png)

**Figure 16. DoubletRemoval GBM131 Naive Violin.**

![Figure 16: DoubletRemoval GBM131 Naive Violin](figures/16_DoubletRemoval_GBM131_Naive_Violin.png)

**Figure 17. DoubletRemoval GBM131 RT2d Scatter.**

![Figure 17: DoubletRemoval GBM131 RT2d Scatter](figures/17_DoubletRemoval_GBM131_RT2d_Scatter.png)

**Figure 18. DoubletRemoval GBM022 Naive Scatter.**

![Figure 18: DoubletRemoval GBM022 Naive Scatter](figures/18_DoubletRemoval_GBM022_Naive_Scatter.png)

**Figure 19. DoubletRemoval GBM022 RT2d Scatter.**

![Figure 19: DoubletRemoval GBM022 RT2d Scatter](figures/19_DoubletRemoval_GBM022_RT2d_Scatter.png)

**Figure 20. DoubletRemoval GBM131 RT2d Violin.**

![Figure 20: DoubletRemoval GBM131 RT2d Violin](figures/20_DoubletRemoval_GBM131_RT2d_Violin.png)

**Figure 21. DoubletRemoval GBM022 RT2d Violin.**

![Figure 21: DoubletRemoval GBM022 RT2d Violin](figures/21_DoubletRemoval_GBM022_RT2d_Violin.png)

**Figure 22. DoubletRemoval RemainingCellCounts.**

![Figure 22: DoubletRemoval RemainingCellCounts](figures/22_DoubletRemoval_RemainingCellCounts.png)

### Highly variable feature selection

Mean–variance plots show 2,000 selected highly variable genes (red) against non-variable genes (black). Several views highlight high-variance genes; the exact sample or object associated with each panel is not recoverable from the screenshots alone.

**Figure 23. HVG Selection Unlabeled View1.**

![Figure 23: HVG Selection Unlabeled View1](figures/23_HVG_Selection_Unlabeled_View1.png)

**Figure 24. HVG Selection Unlabeled View2.**

![Figure 24: HVG Selection Unlabeled View2](figures/24_HVG_Selection_Unlabeled_View2.png)

**Figure 25. HVG Selection Labeled View1.**

![Figure 25: HVG Selection Labeled View1](figures/25_HVG_Selection_Labeled_View1.png)

**Figure 26. HVG Selection Labeled View2.**

![Figure 26: HVG Selection Labeled View2](figures/26_HVG_Selection_Labeled_View2.png)

**Figure 27. HVG Selection Labeled View3.**

![Figure 27: HVG Selection Labeled View3](figures/27_HVG_Selection_Labeled_View3.png)

### PCA and sample-level UMAP exploration

An elbow plot, four UMAP clusterings with different cluster counts, and heatmaps of genes associated with principal components document exploratory dimensionality reduction. The four UMAPs must not be interpreted as a chronological progression or as particular samples without matching R code.

**Figure 28. PCA ElbowPlot.**

![Figure 28: PCA ElbowPlot](figures/28_PCA_ElbowPlot.png)

**Figure 29. UMAP Clustering 9Clusters.**

![Figure 29: UMAP Clustering 9Clusters](figures/29_UMAP_Clustering_9Clusters.png)

**Figure 30. UMAP Clustering 8Clusters.**

![Figure 30: UMAP Clustering 8Clusters](figures/30_UMAP_Clustering_8Clusters.png)

**Figure 31. UMAP Clustering 6Clusters.**

![Figure 31: UMAP Clustering 6Clusters](figures/31_UMAP_Clustering_6Clusters.png)

**Figure 32. UMAP Clustering 4Clusters.**

![Figure 32: UMAP Clustering 4Clusters](figures/32_UMAP_Clustering_4Clusters.png)

**Figure 33. PCA Top10Components Heatmaps Set1.**

![Figure 33: PCA Top10Components Heatmaps Set1](figures/33_PCA_Top10Components_Heatmaps_Set1.png)

**Figure 34. PCA PC1 GeneLoadings Heatmap.**

![Figure 34: PCA PC1 GeneLoadings Heatmap](figures/34_PCA_PC1_GeneLoadings_Heatmap.png)

**Figure 35. PCA Top10Components Heatmaps Set2.**

![Figure 35: PCA Top10Components Heatmaps Set2](figures/35_PCA_Top10Components_Heatmaps_Set2.png)

### Reference-based cell annotation

SingleR-derived labels are projected onto UMAP embeddings. These are computational reference matches, not validated identities: malignant GBM cultures can resemble several normal or developmental reference populations.

**Figure 36. SingleR Annotation UMAP View1.**

![Figure 36: SingleR Annotation UMAP View1](figures/36_SingleR_Annotation_UMAP_View1.png)

**Figure 37. SingleR Annotation UMAP View2.**

![Figure 37: SingleR Annotation UMAP View2](figures/37_SingleR_Annotation_UMAP_View2.png)

**Figure 38. SingleR Annotation UMAP View3.**

![Figure 38: SingleR Annotation UMAP View3](figures/38_SingleR_Annotation_UMAP_View3.png)

**Figure 39. SingleR Annotation UMAP View4.**

![Figure 39: SingleR Annotation UMAP View4](figures/39_SingleR_Annotation_UMAP_View4.png)

### Integrated cell landscape

UMAP embeddings display sample of origin, treatment condition, cluster identities, and SingleR labels. The clear separation by patient-derived culture should be considered when interpreting treatment-associated differences and batch integration.

**Figure 40. IntegratedUMAP BySample.**

![Figure 40: IntegratedUMAP BySample](figures/40_IntegratedUMAP_BySample.png)

**Figure 41. IntegratedUMAP ByCondition.**

![Figure 41: IntegratedUMAP ByCondition](figures/41_IntegratedUMAP_ByCondition.png)

**Figure 42. IntegratedUMAP ByCluster.**

![Figure 42: IntegratedUMAP ByCluster](figures/42_IntegratedUMAP_ByCluster.png)

**Figure 43. IntegratedUMAP SingleR Annotations.**

![Figure 43: IntegratedUMAP SingleR Annotations](figures/43_IntegratedUMAP_SingleR_Annotations.png)

### GBM-associated expression programs

Dot, feature, and violin plots visualize selected GBM markers including SOX2, OLIG2, PDGFRA, EGFR, CD44, CHI3L1, VIM, and NES. Marker expression alone cannot establish definitive GBM cell states.

**Figure 44. GBMStateMarkers Cluster DotPlot.**

![Figure 44: GBMStateMarkers Cluster DotPlot](figures/44_GBMStateMarkers_Cluster_DotPlot.png)

**Figure 45. GBMStateMarkers UMAP FeaturePlots.**

![Figure 45: GBMStateMarkers UMAP FeaturePlots](figures/45_GBMStateMarkers_UMAP_FeaturePlots.png)

**Figure 46. GBMStateMarkers NaiveVsRT Violin.**

![Figure 46: GBMStateMarkers NaiveVsRT Violin](figures/46_GBMStateMarkers_NaiveVsRT_Violin.png)

### Cluster composition and Milo differential abundance

Neighborhood sizes, sample-level cluster proportions, and exploratory Milo differential-abundance statistics are shown. There are only two paired patient-derived cultures, so apparent statistical significance requires careful validation.

**Figure 47. Milo NeighborhoodSize Distribution.**

![Figure 47: Milo NeighborhoodSize Distribution](figures/47_Milo_NeighborhoodSize_Distribution.png)

**Figure 48. ClusterComposition BySample.**

![Figure 48: ClusterComposition BySample](figures/48_ClusterComposition_BySample.png)

**Figure 49. Milo DifferentialAbundance Exploratory.**

![Figure 49: Milo DifferentialAbundance Exploratory](figures/49_Milo_DifferentialAbundance_Exploratory.png)

### Paired pseudobulk differential expression

Pseudobulk library sizes and RT-versus-Naive volcano/MA plots show candidate radiation-associated genes. The labeled plots highlight examples such as GEM and KRT6A; the count chart displays four upregulated and one downregulated genes under the plotted criteria. With two pairs, these are exploratory findings.

**Figure 50. Pseudobulk LibrarySizes.**

![Figure 50: Pseudobulk LibrarySizes](figures/50_Pseudobulk_LibrarySizes.png)

**Figure 51. Pseudobulk RTvsNaive Volcano.**

![Figure 51: Pseudobulk RTvsNaive Volcano](figures/51_Pseudobulk_RTvsNaive_Volcano.png)

**Figure 52. Pseudobulk RTvsNaive MAPlot.**

![Figure 52: Pseudobulk RTvsNaive MAPlot](figures/52_Pseudobulk_RTvsNaive_MAPlot.png)

**Figure 53. Pseudobulk RTvsNaive LabeledVolcano View1.**

![Figure 53: Pseudobulk RTvsNaive LabeledVolcano View1](figures/53_Pseudobulk_RTvsNaive_LabeledVolcano_View1.png)

**Figure 54. Pseudobulk RTvsNaive LabeledVolcano View2.**

![Figure 54: Pseudobulk RTvsNaive LabeledVolcano View2](figures/54_Pseudobulk_RTvsNaive_LabeledVolcano_View2.png)

**Figure 55. Pseudobulk DifferentialGeneCounts.**

![Figure 55: Pseudobulk DifferentialGeneCounts](figures/55_Pseudobulk_DifferentialGeneCounts.png)

### Inferred cluster communication and pathway enrichment

A cluster-to-cluster communication network and GO Biological Process/KEGG enrichment charts are included. Several enrichment categories show gene counts of one, and some labels overlap, limiting interpretability. The network is computationally inferred rather than direct experimental evidence of signaling.

**Figure 56. CellChat ClusterCommunication Network.**

![Figure 56: CellChat ClusterCommunication Network](figures/56_CellChat_ClusterCommunication_Network.png)

**Figure 57. GO BP Upregulated Enrichment.**

![Figure 57: GO BP Upregulated Enrichment](figures/57_GO_BP_Upregulated_Enrichment.png)

**Figure 58. GO BP Downregulated Enrichment.**

![Figure 58: GO BP Downregulated Enrichment](figures/58_GO_BP_Downregulated_Enrichment.png)

**Figure 59. KEGG Upregulated Enrichment.**

![Figure 59: KEGG Upregulated Enrichment](figures/59_KEGG_Upregulated_Enrichment.png)

**Figure 60. KEGG Downregulated Enrichment.**

![Figure 60: KEGG Downregulated Enrichment](figures/60_KEGG_Downregulated_Enrichment.png)

### Mesenchymal-like state and radiation-response genes

Paired MES-like module-score means, a radiation-response marker dot plot, and a scaled-expression heatmap are shown. Both patient-derived cultures have higher displayed MES-like mean scores after RT, but this is descriptive and cannot establish a generalizable effect or causal transition.

**Figure 61. GBM MESlike PairedPatientMeans.**

![Figure 61: GBM MESlike PairedPatientMeans](figures/61_GBM_MESlike_PairedPatientMeans.png)

**Figure 62. RadiationResponse Genes BySample DotPlot.**

![Figure 62: RadiationResponse Genes BySample DotPlot](figures/62_RadiationResponse_Genes_BySample_DotPlot.png)

**Figure 63. RadiationResponse Genes ScaledExpression Heatmap.**

![Figure 63: RadiationResponse Genes ScaledExpression Heatmap](figures/63_RadiationResponse_Genes_ScaledExpression_Heatmap.png)

## Discussion

The analysis is designed to investigate heterogeneity in cultured GBM cells and candidate transcriptional changes following irradiation. Any observed differences must be evaluated against patient-specific variation, cluster composition, sequencing quality, and the very small number of independent samples. Pathway enrichment and ligand–receptor predictions identify hypotheses rather than mechanisms established by experiment.

**The figure-level descriptions above are grounded in the supplied screenshots.** The original numerical tables and complete model diagnostics were not provided, so biological and statistical claims remain provisional.

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

This repository documents an R-based exploratory workflow for examining single-cell transcriptional heterogeneity and radiation-associated candidate responses in two patient-derived glioblastoma cultures. The outputs can guide further validation, but the current dataset is insufficient for definitive claims about radiation resistance mechanisms or clinical biomarkers.

## Repository Structure

```text
GBM-Radiation-scRNAseq/
├── README.md
├── GBM_scRNAseq_Analysis.R
├── requirements.R
├── .gitignore
└── figures/                  # 63 supplied analysis screenshots, embedded in README
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

**Youssef Basem**

**GitHub: @youssefbasem2026-lab**  
GitHub: [@youssefbasem2026-lab](https://github.com/youssefbasem2026-lab)

## License

No license has been assigned yet. A code license can be selected after confirming the intended reuse permissions. The original dataset retains its own attribution and applicable terms.
