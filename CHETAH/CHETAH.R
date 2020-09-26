library(CHETAH)

# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# ref_ndata represents the normolized dgCMatrix object using Seurat from HCL or MCA. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# ref_celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# Create sce object
ref_celltype<- data.frame(celltype = ref_celltype$celltype,stringsAsFactors = F)
rownames(ref_celltype)<- colnames(ref_ndata)
ref_sce <- SingleCellExperiment(assays = list(counts = as.matrix(ref_ndata)),
                                colData = ref_celltype)

test_sce <- SingleCellExperiment(assays = list(counts = as.matrix(test_data)))

# predict
test_sce <- CHETAHclassifier(input = test_sce,
                             ref_cells = ref_sce)

celltype$chetah<- test_sce$celltype_CHETAH