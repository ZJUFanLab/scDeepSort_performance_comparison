library(singleCellNet)
# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# ref_ndata represents the normolized dgCMatrix object using Seurat from HCL or MCA. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# ref_celltype represents the data.frame object containg two columns, namely cell barcode and cell type.


colnames(ref_celltype)<- c('cell','ann')
colnames(celltype)<- c('cell','ann')
ref_ndata<- ref_ndata[rownames(ref_ndata) %in% rownames(ndata),]
#train
class_info<- scn_train(stTrain = ref_celltype,
                      expTrain = ref_ndata,
                      dLevel = "ann",
                      colName_samp = "cell")
# Predict
crParkall<- scn_predict(class_info[['cnProc']], ndata, nrand = 2)

stPark <- get_cate(classRes = crParkall,
                   sampTab = celltype, dLevel = "ann", sid = "cell",
                   nrand = 50)