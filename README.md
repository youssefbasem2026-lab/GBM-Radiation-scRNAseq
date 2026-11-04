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

<img width="907" height="342" alt="Screenshot 2026-11-05 003048" src="https://github.com/user-attachments/assets/8f4ef90a-d7a7-4e58-a511-0f29294d747b" />

**Figure 3. QC PreFilter GBM022 Naive.**

<img width="911" height="327" alt="Screenshot 2026-11-05 003054" src="https://github.com/user-attachments/assets/358739b7-2d38-483a-b98d-5f4f7140ba1a" />

**Figure 4. QC PreFilter GBM022 RT2d.**

<img width="887" height="292" alt="Screenshot 2026-11-05 003059" src="https://github.com/user-attachments/assets/f49214a7-42b2-4a60-aae0-5485157351b6" />

**Figure 5. QC PostFilter GBM131 Naive Overview.**

<img width="876" height="267" alt="Screenshot 2026-11-05 003105" src="https://github.com/user-attachments/assets/39ff1570-3997-4dbf-9f97-e1b615f493e6" />

**Figure 6. QC PostFilter GBM131 Naive Refined.**

<img width="841" height="323" alt="Screenshot 2026-11-05 003112" src="https://github.com/user-attachments/assets/accc6ac4-2cc9-4714-8e5d-7745b6721cd6" />

**Figure 7. QC PostFilter GBM131 RT2d.**

<img width="817" height="267" alt="Screenshot 2026-11-05 003119" src="https://github.com/user-attachments/assets/eb68a678-4c6f-400d-8ce8-275dd443595f" />

**Figure 8. QC PostFilter GBM022 Naive.**

<img width="823" height="336" alt="Screenshot 2026-11-05 003125" src="https://github.com/user-attachments/assets/05914b68-d7fa-421e-9bbf-fac35b9af6ac" />

**Figure 9. QC PostFilter GBM022 RT2d.**

<img width="810" height="306" alt="Screenshot 2026-11-05 003134" src="https://github.com/user-attachments/assets/71e85f77-ab7d-4982-a313-18f3df3a5b3d" />

### Doublet detection and removal

scDblFinder doublet scores and class assignments are visualized before and after removal. The retained-cell bar chart compares the four samples; bar heights are descriptive and are not a formal treatment-effect test.

**Figure 10. DoubletDetection GBM131 Naive Scatter.**

<img width="782" height="267" alt="Screenshot 2026-11-05 003140" src="https://github.com/user-attachments/assets/b8563e37-9da4-47c3-9bb0-2c8be4ce0b52" />

**Figure 11. DoubletDetection GBM131 Naive Violin.**

<img width="786" height="286" alt="Screenshot 2026-11-05 003145" src="https://github.com/user-attachments/assets/943ecea3-ae92-464e-bed7-d55861465102" />

**Figure 12. DoubletDetection GBM131 RT2d Scatter.**

<img width="767" height="298" alt="Screenshot 2026-11-05 003150" src="https://github.com/user-attachments/assets/cfd21fdc-e6dc-4ab2-aa6f-2f0d41808f17" />

**Figure 13. DoubletDetection GBM131 RT2d Violin.**

<img width="751" height="267" alt="Screenshot 2026-11-05 003155" src="https://github.com/user-attachments/assets/b01f5683-d70d-4508-81be-d2f1fd7edff2" />

**Figure 14. DoubletDetection GBM022 Naive Scatter.**

<img width="742" height="282" alt="Screenshot 2026-11-05 003201" src="https://github.com/user-attachments/assets/c50a45e0-bbd6-43df-be86-1b84b44ca263" />

**Figure 15. DoubletRemoval GBM131 Naive Scatter.**

<img width="738" height="282" alt="Screenshot 2026-11-05 003206" src="https://github.com/user-attachments/assets/0531a013-4196-4a53-870d-c806df7e71b5" />

**Figure 16. DoubletRemoval GBM131 Naive Violin.**

<img width="718" height="280" alt="Screenshot 2026-11-05 003212" src="https://github.com/user-attachments/assets/389dc25e-ff14-47f8-a03c-9d081ff30016" />

**Figure 17. DoubletRemoval GBM131 RT2d Scatter.**

<img width="695" height="297" alt="Screenshot 2026-11-05 003218" src="https://github.com/user-attachments/assets/67c7b5b2-0e1a-45c7-bd48-952dcaa168eb" />

**Figure 18. DoubletRemoval GBM022 Naive Scatter.**

<img width="686" height="272" alt="Screenshot 2026-11-05 003225" src="https://github.com/user-attachments/assets/dd1ebc6c-5ac6-450f-99e6-34115b31ea93" />

**Figure 19. DoubletRemoval GBM022 RT2d Scatter.**

<img width="672" height="305" alt="Screenshot 2026-11-05 003230" src="https://github.com/user-attachments/assets/7ecab579-3eca-4ba7-9a8a-1e64a79465f0" />

**Figure 20. DoubletRemoval GBM131 RT2d Violin.**

<img width="666" height="277" alt="Screenshot 2026-11-05 003238" src="https://github.com/user-attachments/assets/2ff82852-7796-4d44-a904-2f5ec65366d1" />

**Figure 21. DoubletRemoval GBM022 RT2d Violin.**

<img width="652" height="267" alt="Screenshot 2026-11-05 003243" src="https://github.com/user-attachments/assets/d83926da-9aec-43b7-bf47-90bd15de2767" />

**Figure 22. DoubletRemoval RemainingCellCounts.**

<img width="645" height="297" alt="Screenshot 2026-11-05 003249" src="https://github.com/user-attachments/assets/938ed0f8-a385-4837-8e4f-7a53760301c4" />

### Highly variable feature selection

Mean–variance plots show 2,000 selected highly variable genes (red) against non-variable genes (black). Several views highlight high-variance genes; the exact sample or object associated with each panel is not recoverable from the screenshots alone.

**Figure 23. HVG Selection Unlabeled View1.**

<img width="627" height="305" alt="Screenshot 2026-11-05 003254" src="https://github.com/user-attachments/assets/7dbfba88-521a-4f6f-9f39-bb93305dfaa5" />

**Figure 24. HVG Selection Unlabeled View2.**

<img width="622" height="275" alt="Screenshot 2026-11-05 003258" src="https://github.com/user-attachments/assets/5b086a4c-d1af-4e6f-946e-4ea09a299a26" />

**Figure 25. HVG Selection Labeled View1.**

<img width="612" height="278" alt="Screenshot 2026-11-05 003303" src="https://github.com/user-attachments/assets/bf81bf66-5fe9-4b07-ada3-1004c7f5926c" />

**Figure 26. HVG Selection Labeled View2.**

<img width="597" height="326" alt="Screenshot 2026-11-05 003309" src="https://github.com/user-attachments/assets/997d7658-ac0e-48d1-966f-f7637967baa6" />

**Figure 27. HVG Selection Labeled View3.**

<img width="675" height="321" alt="Screenshot 2026-11-05 003314" src="https://github.com/user-attachments/assets/a43df47c-c24b-465d-8540-a0214379ac62" />

### PCA and sample-level UMAP exploration

An elbow plot, four UMAP clusterings with different cluster counts, and heatmaps of genes associated with principal components document exploratory dimensionality reduction. The four UMAPs must not be interpreted as a chronological progression or as particular samples without matching R code.

**Figure 28. PCA ElbowPlot.**

<img width="660" height="291" alt="Screenshot 2026-11-05 003320" src="https://github.com/user-attachments/assets/06310ab2-8b97-4742-94d5-04ad1216cc11" />

**Figure 29. UMAP Clustering 9Clusters.**

<img width="647" height="311" alt="Screenshot 2026-11-05 003325" src="https://github.com/user-attachments/assets/319ab3f7-c4c1-4f8f-bfa9-fb7398a832c1" />

**Figure 30. UMAP Clustering 8Clusters.**

<img width="640" height="318" alt="Screenshot 2026-11-05 003330" src="https://github.com/user-attachments/assets/e69e0fb1-ae0b-4fc5-a763-8fc13deb6914" />

**Figure 31. UMAP Clustering 6Clusters.**

<img width="627" height="342" alt="Screenshot 2026-11-05 003336" src="https://github.com/user-attachments/assets/751d505d-f39a-422a-b10e-94adbca116c1" />

**Figure 32. UMAP Clustering 4Clusters.**

<img width="617" height="338" alt="Screenshot 2026-11-05 003346" src="https://github.com/user-attachments/assets/b88a1ec1-2f4b-4770-8601-7de6272b8668" />

**Figure 33. PCA Top10Components Heatmaps Set1.**

<img width="738" height="372" alt="Screenshot 2026-11-05 003350" src="https://github.com/user-attachments/assets/69e228d4-03ca-41fd-b14c-eb8d82cfc2b8" />

**Figure 34. PCA PC1 GeneLoadings Heatmap.**

<img width="247" height="97" alt="Screenshot 2026-11-05 003355" src="https://github.com/user-attachments/assets/b20c2342-e4a7-4717-bb4a-a926e1e650be" />

**Figure 35. PCA Top10Components Heatmaps Set2.**

<img width="732" height="346" alt="Screenshot 2026-11-05 003403" src="https://github.com/user-attachments/assets/238a0a50-6c87-440c-a8ad-f306c839ad79" />

### Reference-based cell annotation

SingleR-derived labels are projected onto UMAP embeddings. These are computational reference matches, not validated identities: malignant GBM cultures can resemble several normal or developmental reference populations.

**Figure 36. SingleR Annotation UMAP View1.**

<img width="705" height="331" alt="Screenshot 2026-11-05 003408" src="https://github.com/user-attachments/assets/efad1884-91fd-4082-abd8-faa388d30182" />

**Figure 37. SingleR Annotation UMAP View2.**

<img width="700" height="327" alt="Screenshot 2026-11-05 003412" src="https://github.com/user-attachments/assets/ba27c0cf-c573-42f9-9e07-4cd6fb0fb380" />

**Figure 38. SingleR Annotation UMAP View3.**

<img width="688" height="347" alt="Screenshot 2026-11-05 003417" src="https://github.com/user-attachments/assets/2e7fe9ed-47bf-4fdb-a8f2-0d89aae270a5" />

**Figure 39. SingleR Annotation UMAP View4.**

<img width="676" height="321" alt="Screenshot 2026-11-05 003422" src="https://github.com/user-attachments/assets/7f5029f4-fb36-44c4-87ce-227b9d4d9fe2" />

### Integrated cell landscape

UMAP embeddings display sample of origin, treatment condition, cluster identities, and SingleR labels. The clear separation by patient-derived culture should be considered when interpreting treatment-associated differences and batch integration.

**Figure 40. IntegratedUMAP BySample.**

<img width="668" height="300" alt="Screenshot 2026-11-05 003427" src="https://github.com/user-attachments/assets/fcb6017e-83e3-421f-9dca-1fd57ba3b560" />

**Figure 41. IntegratedUMAP ByCondition.**

<img width="651" height="267" alt="Screenshot 2026-11-05 003433" src="https://github.com/user-attachments/assets/412586d7-d05f-4890-8bda-d9eb764ed013" />

**Figure 42. IntegratedUMAP ByCluster.**

<img width="635" height="252" alt="Screenshot 2026-11-05 003438" src="https://github.com/user-attachments/assets/a7fc64a2-af74-4d12-9553-e1dc8c3ee56a" />

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
