library(SingleCellExperiment)
library(scmap)

# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# ref_ndata represents the normolized dgCMatrix object using Seurat from HCL or MCA. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# ref_celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# Create ref_ndata sce object
ref_celltype<- data.frame(cell_type1 = ref_celltype$celltype,stringsAsFactors = F)
rownames(ref_celltype)<- colnames(ref_ndata)
ref_ndata_sce<- SingleCellExperiment(assays = list(normcounts = as.matrix(ref_ndata)),colData = ref_celltype)
logcounts(ref_ndata_sce) <- normcounts(ref_ndata_sce)
rowData(ref_ndata_sce)$feature_symbol <- rownames(ref_ndata_sce)
isSpike(ref_ndata_sce, "ERCC") <- grepl("^ERCC-", rownames(ref_ndata_sce))
ref_ndata_sce <- ref_ndata_sce[!duplicated(rownames(ref_ndata_sce)), ]
ref_ndata_sce <- selectFeatures(ref_ndata_sce, suppress_plot = F)

# Create test data sce object
test_ndata_sce<- SingleCellExperiment(assays = list(normcounts = as.matrix(ndata)))
logcounts(test_ndata_sce) <- normcounts(test_ndata_sce)
rowData(test_ndata_sce)$feature_symbol <- rownames(test_ndata_sce)
isSpike(test_ndata_sce, "ERCC") <- grepl("^ERCC-", rownames(test_ndata_sce))
test_ndata_sce <- test_ndata_sce[!duplicated(rownames(test_ndata_sce)), ]
test_ndata_sce <- selectFeatures(test_ndata_sce, suppress_plot = FALSE)

# Annotating cell-based
ref_ndata_sce <- indexCell(ref_ndata_sce)
scmapCell_results <- scmapCell(
  test_ndata_sce, 
  list(
    scmap = metadata(ref_ndata_sce)$scmap_cell_index
  )
)
scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(ref_ndata_sce)$cell_type1)
  )
)
celltype$scmap_cell<- scmapCell_clusters$scmap_cluster_labs


# Annotating cluster-based
ref_ndata_sce <- indexCluster(ref_ndata_sce)
scmapCluster_results <- scmapCluster(
  test_ndata_sce, 
  list(
    scmap = metadata(ref_ndata_sce)$scmap_cluster_index
  )
)
celltype$scmap_cluster<- scmapCluster_results$scmap_cluster_labs
