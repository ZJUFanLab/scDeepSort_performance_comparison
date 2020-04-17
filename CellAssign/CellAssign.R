library(SingleCellExperiment)
library(cellassign)
library(scran)

# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.
# marker_file is obtained from [marker_file_making.R]
marker_file<- marker_file[rownames(marker_file) %in% rownames(ndata),]
cellname<- rep(1,ncol(marker_file))
for (i in 1:ncol(marker_file)) {
  d1<- marker_file[,i]
  if (all(d1 == 0)) {
    cellname[i]<- 0
  }
}
cellname<- which(cellname == 1)
marker_file<- marker_file[,cellname]

# making colData and rowData
ndata_colData<- data.frame(Cell = celltype$cell_barcode,stringsAsFactors = F)
rownames(ndata_colData)<- ndata_colData$Cell

ndata_rowData <- data.frame(Gene = rownames(ndata),stringsAsFactors = F)
rownames(ndata_rowData) <- ndata_rowData$Gene

# Create sce object
ndata_sce <- SingleCellExperiment(assays = list(counts = as.matrix(ndata)),colData = ndata_colData,rowData = ndata_rowData)
ndata_sce <- computeSumFactors(ndata_sce)
ndata_sce_Size_factors <- sizeFactors(ndata_sce)

# cellassign
ndata_cellassign_fit <- cellassign(exprs_obj = ndata_sce[rownames(marker_file),], 
                                   marker_gene_info = as.matrix(marker_file), 
                                   s = ndata_sce_Size_factors, 
                                   learning_rate = 1e-2, 
                                   shrinkage = TRUE,
                                   verbose = TRUE)
