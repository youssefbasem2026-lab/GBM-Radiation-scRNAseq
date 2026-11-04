# Packages used in the analysis script.
# Run installation commands selectively; check R/Bioconductor compatibility first.
cran_packages <- c('Seurat', 'ggplot2', 'dplyr', 'patchwork', 'harmony', 'ggrepel', 'scales')
bioconductor_packages <- c('SingleCellExperiment', 'scDblFinder', 'SingleR', 'celldex',
                          'miloR', 'edgeR', 'clusterProfiler', 'org.Hs.eg.db', 'enrichplot')
optional_packages <- c('CellChat', 'slingshot')

# install.packages(cran_packages)
# install.packages('BiocManager')
# BiocManager::install(bioconductor_packages, ask = FALSE, update = FALSE)
# CellChat installation may require a separate installation method/version.
# BiocManager::install('slingshot', ask = FALSE, update = FALSE)
