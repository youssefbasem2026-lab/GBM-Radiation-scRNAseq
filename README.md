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

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/3d13ee6a-37c4-4d2d-a6ef-0b06e167926d" />


**Figure 2. QC PreFilter GBM131 RT2d.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/e91fa93e-888f-45fb-b17f-87399607122a" />

**Figure 3. QC PreFilter GBM022 Naive.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/bad8e8ab-c450-4078-a136-f3b1e0e304ff" />

**Figure 4. QC PreFilter GBM022 RT2d.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/ee979896-6c11-4b42-a17b-66f40eef55d3" />

**Figure 5. QC PostFilter GBM131 Naive Overview.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/dbfdae6c-ebce-4b2c-824b-d5fc587db1c4" />

**Figure 6. QC PostFilter GBM131 Naive Refined.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/dbe86305-97e3-48c2-8e4a-aa1ced550e74" />

**Figure 7. QC PostFilter GBM131 RT2d.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/1e134c55-4f18-4db3-afe6-b682567af298" />

**Figure 8. QC PostFilter GBM022 Naive.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/bbff3e84-50df-43c5-a9cc-3022c2ee7a27" />

**Figure 9. QC PostFilter GBM022 RT2d.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/4ce10d80-1a9d-4f5c-8dfb-8b10dc706ec5" />

### Doublet detection and removal

scDblFinder doublet scores and class assignments are visualized before and after removal. The retained-cell bar chart compares the four samples; bar heights are descriptive and are not a formal treatment-effect test.

**Figure 10. DoubletDetection GBM131 Naive Scatter.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/caaae69a-af53-4312-9b81-44bc8750340d" />

**Figure 11. DoubletDetection GBM131 Naive Violin.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/567d08d1-ec8e-4e74-874e-824db4fa3ccb" />

**Figure 12. DoubletDetection GBM131 RT2d Scatter.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/28818d42-84f6-4e04-8ada-3b07db0a557c" />

**Figure 13. DoubletDetection GBM131 RT2d Violin.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/af4b71bc-e80f-4401-ba73-886271cbe813" />

**Figure 14. DoubletDetection GBM022 Naive Scatter.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/95c07352-af43-48e9-bd72-0322f5ace664" />

**Figure 15. DoubletRemoval GBM131 Naive Scatter.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/ce79262b-5d86-4485-bae4-0ee30590c0e2" />

**Figure 16. DoubletRemoval GBM131 RT2d Scatter.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/eead246d-3861-4d4d-98fa-b37e2fae7f66" />

**Figure 17. DoubletRemoval GBM022 Naive Scatter.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/4a35865e-3bfc-4f5a-996f-4c6bd4e74be2" />

**Figure 18. DoubletRemoval GBM022 RT2d Scatter.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/f57b70f0-9b14-424f-acc3-155e901de6c1" />

**Figure 19. DoubletRemoval GBM131 RT2d Violin.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/cabf49ed-33f1-4c62-96ab-19e437cb3b18" />

**Figure 20. DoubletRemoval GBM022 RT2d Violin.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/03a1e08e-a0f6-486c-8de8-ecbbbda30d34" />

**Figure 21. DoubletRemoval RemainingCellCounts.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/9e31c1e4-4648-4b76-bccb-36361597f534" />

### Highly variable feature selection

Mean–variance plots show 2,000 selected highly variable genes (red) against non-variable genes (black). Several views highlight high-variance genes; the exact sample or object associated with each panel is not recoverable from the screenshots alone.

**Figure 22. HVG Selection Unlabeled View1.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/98680ba9-f7b9-4348-a88f-cca0b73f5e02" />

**Figure 23. HVG Selection Unlabeled View2.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/706aaa90-c706-4a86-9f67-8583d199c948" />

**Figure 24. HVG Selection Labeled View1.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/ba9df3d2-b57c-4f90-9320-be65be647b75" />

**Figure 25. HVG Selection Labeled View2.**
<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/97b1ac66-9c64-4c0a-9f09-50cb36097f07" />


**Figure 26,27. HVG Selection Labeled View3.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/aeabfa1b-97b8-4e41-b736-c8e3ef79b5b2" />
<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/546f1d91-4df5-46c6-bfa5-052495cae0b8" />

### PCA and sample-level UMAP exploration

An elbow plot, four UMAP clusterings with different cluster counts, and heatmaps of genes associated with principal components document exploratory dimensionality reduction. The four UMAPs must not be interpreted as a chronological progression or as particular samples without matching R code.

**Figure 28. PCA ElbowPlot.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/edb08b7e-2e5b-4685-bdce-074e84e5a4ae" />

**Figure 29. UMAP Clustering 9Clusters.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/03e91c7b-9c13-41c3-9b14-d7e0c26f952e" />

**Figure 30. UMAP Clustering 8Clusters.**
<img width="640" height="318" alt="Screenshot 2026-11-05 003330" src="https://github.com/user-attachments/assets/07898542-6653-4fbe-9b11-4238896ecb86" />

**Figure 31. UMAP Clustering 6Clusters.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/5240e9aa-48b5-4afe-8320-08d15bd685e6" />

**Figure 32. UMAP Clustering 4Clusters.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/13d7aa47-03f2-4361-9f8d-2fd397bfee99" />

**Figure 33. PCA Top10Components Heatmaps Set1.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/1a7725e8-730e-45f7-9d23-094c6c7c809b" />
<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/201c299c-428f-4ecc-912e-92be7645f051" />

**Figure 34. PCA PC1 GeneLoadings Heatmap.**

<img width="247" height="97" alt="Screenshot 2026-11-05 003355" src="https://github.com/user-attachments/assets/812aadf6-4697-45d0-93a8-144b53d365ee" />

**Figure 35. PCA Top10Components Heatmaps Set2.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/e4329008-e97f-4d45-894e-b09f93e33cdc" />

### Reference-based cell annotation

SingleR-derived labels are projected onto UMAP embeddings. These are computational reference matches, not validated identities: malignant GBM cultures can resemble several normal or developmental reference populations.

**Figure 36. SingleR Annotation UMAP View1.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/12089310-e040-4e7a-808c-93e6dc7a910b" />

**Figure 37. SingleR Annotation UMAP View2.**

<img width="700" height="327" alt="Screenshot 2026-11-05 003412" src="https://github.com/user-attachments/assets/9dcea30d-b987-4e8e-9044-49a99ce2ef53" />

**Figure 38. SingleR Annotation UMAP View3.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/80a37112-93bf-4565-b3c6-27dadaf1ff3c" />

**Figure 39. SingleR Annotation UMAP View4.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/83d58cea-ba64-405a-bcbd-ecba8683b6ca" />

### Integrated cell landscape

UMAP embeddings display sample of origin, treatment condition, cluster identities, and SingleR labels. The clear separation by patient-derived culture should be considered when interpreting treatment-associated differences and batch integration.

**Figure 40. IntegratedUMAP BySample.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/49e9322e-b2d6-4ccd-8442-002754952f89" />

**Figure 41. IntegratedUMAP ByCondition.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/75f9285d-8b14-4d28-af81-e35695c34698" />

**Figure 42. IntegratedUMAP ByCluster.**

<img width="927" height="483" alt="image" src="https://github.com/user-attachments/assets/2d372878-2d8a-4dd5-a9d4-3f4edbd805d6" />

**Figure 43. IntegratedUMAP SingleR Annotations.**

<img width="631" height="292" alt="Screenshot 2026-11-05 003443" src="https://github.com/user-attachments/assets/c9d15834-f223-47a9-a04c-671ee20255d2" />

### GBM-associated expression programs

Dot, feature, and violin plots visualize selected GBM markers including SOX2, OLIG2, PDGFRA, EGFR, CD44, CHI3L1, VIM, and NES. Marker expression alone cannot establish definitive GBM cell states.

**Figure 44. GBMStateMarkers Cluster DotPlot.**

<img width="652" height="318" alt="Screenshot 2026-11-05 003449" src="https://github.com/user-attachments/assets/5b8f12f3-b794-4e79-90fa-d79fee8cbbf1" />

**Figure 45. GBMStateMarkers UMAP FeaturePlots.**

<img width="687" height="325" alt="Screenshot 2026-11-05 003454" src="https://github.com/user-attachments/assets/700f272b-21b4-47e0-ba40-5a0192c4c2d5" />

**Figure 46. GBMStateMarkers NaiveVsRT Violin.**

<img width="672" height="291" alt="Screenshot 2026-11-05 003458" src="https://github.com/user-attachments/assets/d15f3c75-c271-44b5-ae01-ac67f63c2436" />

### Cluster composition and Milo differential abundance

Neighborhood sizes, sample-level cluster proportions, and exploratory Milo differential-abundance statistics are shown. There are only two paired patient-derived cultures, so apparent statistical significance requires careful validation.

**Figure 47. Milo NeighborhoodSize Distribution.**

<img width="660" height="257" alt="Screenshot 2026-11-05 003504" src="https://github.com/user-attachments/assets/919138ad-c720-4baf-a439-b32c8f1a58aa" />

**Figure 48. ClusterComposition BySample.**

<img width="655" height="306" alt="Screenshot 2026-11-05 003508" src="https://github.com/user-attachments/assets/78cbd0bd-a254-44c3-82f5-6566b76d962c" />

**Figure 49. Milo DifferentialAbundance Exploratory.**

<img width="642" height="280" alt="Screenshot 2026-11-05 003514" src="https://github.com/user-attachments/assets/ffb94f01-23e3-4d95-a989-2849de318a01" />

### Paired pseudobulk differential expression

Pseudobulk library sizes and RT-versus-Naive volcano/MA plots show candidate radiation-associated genes. The labeled plots highlight examples such as GEM and KRT6A; the count chart displays four upregulated and one downregulated genes under the plotted criteria. With two pairs, these are exploratory findings.

**Figure 50. Pseudobulk LibrarySizes.**

<img width="637" height="297" alt="Screenshot 2026-11-05 003521" src="https://github.com/user-attachments/assets/3b89a218-e88d-4dbe-81e9-5888d04f8c32" />

**Figure 51. Pseudobulk RTvsNaive Volcano.**

<img width="628" height="272" alt="Screenshot 2026-11-05 003526" src="https://github.com/user-attachments/assets/04fd4845-1ae6-4231-932b-e258fa4a8b5f" />

**Figure 52. Pseudobulk RTvsNaive MAPlot.**

<img width="610" height="300" alt="Screenshot 2026-11-05 003530" src="https://github.com/user-attachments/assets/2c3714cf-9cfd-4d26-b811-3a48946f3ee5" />

**Figure 53. Pseudobulk RTvsNaive LabeledVolcano View1.**

<img width="650" height="293" alt="Screenshot 2026-11-05 003536" src="https://github.com/user-attachments/assets/b42078d6-55e9-437c-980e-8719524a70dc" />

**Figure 54. Pseudobulk RTvsNaive LabeledVolcano View2.**

<img width="640" height="268" alt="Screenshot 2026-11-05 003542" src="https://github.com/user-attachments/assets/1fb9b34e-3751-4622-a6ef-51aa2ea4cdd1" />

**Figure 55. Pseudobulk DifferentialGeneCounts.**

<img width="631" height="240" alt="Screenshot 2026-11-05 003548" src="https://github.com/user-attachments/assets/042d3bab-2388-4600-ae9c-e1c2e98ac6d7" />

### Inferred cluster communication and pathway enrichment

A cluster-to-cluster communication network and GO Biological Process/KEGG enrichment charts are included. Several enrichment categories show gene counts of one, and some labels overlap, limiting interpretability. The network is computationally inferred rather than direct experimental evidence of signaling.

**Figure 56. CellChat ClusterCommunication Network.**

<img width="797" height="375" alt="Screenshot 2026-11-05 003553" src="https://github.com/user-attachments/assets/ec7699f6-17f6-45b2-bc16-8a5dceacc670" />

**Figure 57. GO BP Upregulated Enrichment.**

<img width="837" height="362" alt="Screenshot 2026-11-05 003558" src="https://github.com/user-attachments/assets/d3a8ef35-523a-4182-8c1a-1f8929dd9bc8" />

**Figure 58. GO BP Downregulated Enrichment.**

<img width="902" height="457" alt="Screenshot 2026-11-05 003602" src="https://github.com/user-attachments/assets/9775dfd2-e100-4cb7-92b9-283262afe337" />

**Figure 59. KEGG Upregulated Enrichment.**

<img width="895" height="428" alt="Screenshot 2026-11-05 003606" src="https://github.com/user-attachments/assets/fe5bd992-e83c-40b2-aaf2-81eb3f3e1616" />

**Figure 60. KEGG Downregulated Enrichment.**

<img width="678" height="320" alt="Screenshot 2026-11-05 003610" src="https://github.com/user-attachments/assets/d13c28d8-f307-42b7-a583-bcda5dfe55d1" />

### Mesenchymal-like state and radiation-response genes

Paired MES-like module-score means, a radiation-response marker dot plot, and a scaled-expression heatmap are shown. Both patient-derived cultures have higher displayed MES-like mean scores after RT, but this is descriptive and cannot establish a generalizable effect or causal transition.

**Figure 61. GBM MESlike PairedPatientMeans.**

<img width="665" height="288" alt="Screenshot 2026-11-05 003615" src="https://github.com/user-attachments/assets/a4feab38-8c5d-4fc6-8d3b-1adc67b20d09" />

**Figure 62. RadiationResponse Genes BySample DotPlot.**

<img width="653" height="282" alt="Screenshot 2026-11-05 003620" src="https://github.com/user-attachments/assets/9ccf3543-211c-45a7-afc0-3508b4a8313c" />

**Figure 63. RadiationResponse Genes ScaledExpression Heatmap.**

<img width="802" height="450" alt="Screenshot 2026-11-05 003624" src="https://github.com/user-attachments/assets/9acbd93b-f113-42c9-a921-cbdb3e0537a4" />

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
