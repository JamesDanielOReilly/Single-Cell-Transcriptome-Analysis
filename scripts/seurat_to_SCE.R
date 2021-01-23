library(scater)
library(Seurat)
library(slalom)

# Convert Seurat object to SingleCellExperiment object
pbmc.rds <- readRDS(file = "/home/james/Documents/leuven/year-2/AMSA/Single-Cell-Transcriptome-Analysis/output/seurat_object.rds")
pbmc.sce <- as.SingleCellExperiment(pbmc.rds)

# Format data correctly for slalom 
logcounts.matrix <- as.matrix(SingleCellExperiment::logcounts(pbmc.sce))
pbmc <- SingleCellExperiment::SingleCellExperiment(
  assays = list(logcounts = logcounts.matrix)
)

# Load geneset
gmtfile <- "/home/james/Documents/leuven/year-2/AMSA/Single-Cell-Transcriptome-Analysis/data/genesets/MSigDB_hallmark.gmt"
genesets <- GSEABase::getGmt(gmtfile)

# Creating the model with a set number of hidden factors and minimum number of genes per gene set
model <- newSlalomModel(pbmc, genesets, n_hidden = 5, min_genes = 10)

# Initialising the model and set seed for reproducible analysis
model <- initSlalom(model, seed = 99)

# Training the model
trained.model <- trainSlalom(model, minIterations = 400, nIterations = 2000, shuffle = TRUE, seed = 99)

# View most relevant terms and their respective gene set sizes
topTerms(model)
plotRelevance(model, mad_filter = 0.1, unannotated_dense = TRUE, unannotated_sparse = FALSE)

# View most relevant annotated and unannotated factors
plotTerms(model, mad_filter = 0.1, unannotated_dense = FALSE, unannotated_sparse = FALSE)
plotTerms(model, mad_filter = 0.1, unannotated_dense = TRUE, unannotated_sparse = FALSE)

# View loadings for specific factors
plotLoadings(model, "hidden01")
plotLoadings(model, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
