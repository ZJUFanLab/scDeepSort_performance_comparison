library(scID)

# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# ref_ndata represents the normolized dgCMatrix object using Seurat from HCL or MCA. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# ref_celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

ref_cluster<- unique(ref_celltype$celltype)
ref_cluster<- data.frame(cluster = 1:length(ref_cluster),celltype = ref_cluster,stringsAsFactors = F)
ref_celltype$cluster<- 'NO'
for (i in 1:nrow(ref_cluster)) {
  ref_celltype[ref_celltype$celltype == ref_cluster$celltype[i],]$cluster<- ref_cluster$cluster[i]
}

ref_Seurat<- CreateSeuratObject(counts = ref_ndata)
Idents(ref_Seurat) <- ref_celltype$cluster
ref_markers<- FindAllMarkers(object = ref_Seurat,test.use = 'MAST',only.pos = F,logfc.threshold = 0.5)
# train and predict
scID_output <- scid_multiclass(target_gem = as.matrix(ndata),
                               reference_gem = as.matrix(ref_ndata),
                               reference_clusters = Idents(ref_Seurat),
                               markers = ref_markers)

celltype$scID<- scID_output$labels


