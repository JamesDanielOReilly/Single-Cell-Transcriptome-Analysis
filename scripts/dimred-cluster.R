# loading libraries
library(dplyr)
library(Seurat)
library(patchwork)

# Load the data
pbmc.data <- Read10X(data.dir = "/home/james/Documents/leuven/year-2/AMSA/Single-Cell-Transcriptome-Analysis/data/10X-genomics/PBMC_3K")

# Initialize the Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Quality Control and cell selection: filtering low-quality cells
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalisation
pbmc <- NormalizeData(pbmc)

# Identification of highly-variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Data Scaling
all.genes.pbmc <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes.pbmc)

# Dimensionality Reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examining the PCA results
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca", combine = TRUE, col = "darkblue")
DimPlot(pbmc, reduction = "pca", label = FALSE, cols = c("darkblue"))

DimHeatmap(pbmc, dims = 1:2, cells = NULL, balanced = TRUE, fast = FALSE, slot = "scale.data")
DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE, fast = FALSE, slot = "scale.data")

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

# Save Seurat object
saveRDS(pbmc, file = "/home/james/Documents/leuven/year-2/AMSA/Single-Cell-Transcriptome-Analysis/output/pbmc.rds")
