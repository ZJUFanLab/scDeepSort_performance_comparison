library(monocle)
library(BiocGenerics)

# For each external testing dataset, we transformed it into a cds object for running Garnett.
# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# processing pdata
pdata<- data.frame(tsne_1 = 0, tsne_2 = 0, Size_Factor = 1, FACS_type = celltype$cell_type,stringsAsFactors = F)
rownames(pdata)<- celltype$cell_barcode
pdata$tsne_1<- sample(x = 0:nrow(pdata),size = nrow(pdata))/nrow(pdata)
pdata$tsne_2<- sample(x = 0:nrow(pdata),size = nrow(pdata))/nrow(pdata)

# processing fdata
fdata<- data.frame(gene_short_name = rownames(ndata),num_cells_expressed = 0,stringsAsFactors = F)
genename<- fdata$gene_short_name
expressed_cells<- NULL
for (j in 1:length(genename)) {
  d1<- as.numeric(ndata[j,])
  expressed_cells[j]<- length(d1[d1 > 0])
}
fdata$num_cells_expressed<- expressed_cells
rownames(fdata)<- fdata$gene_short_name

# making cds object
fdata <- new("AnnotatedDataFrame", data = fdata)
pdata <- new("AnnotatedDataFrame", data = pdata)
ndata_cds <- newCellDataSet(ndata,phenoData = pdata,featureData = fdata)
ndata_cds<- estimateSizeFactors(ndata_cds)