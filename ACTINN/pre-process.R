# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# ref_ndata represents the normolized dgCMatrix object using Seurat from HCL or MCA. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# ref_celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# training data
ref_ndata<- as.matrix(ref_ndata)
write.csv(ndata,file = 'ref_ndata_train.csv')
R.utils::gzip('ref_ndata_train.csv',remove = F)
write.table(ref_celltype,file = 'ref_celltype.txt',
            quote = F,sep = '\t',row.names = F,col.names = F)
R.utils::gzip('ref_celltype.txt',remove = F)

# test data
ndata<- as.matrix(ndata)
write.csv(ndata,file = 'ndata.csv')
R.utils::gzip('ndata.csv',remove = F)
