library(SingleR)
library(Biobase)

# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# ref_ndata represents the normolized dgCMatrix object using Seurat from HCL or MCA. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# ref_celltype represents the data.frame object containg two columns, namely cell barcode and cell type.
# ref_name ('HCL' or 'MCA')
# species ('Human' or 'Mouse')

# Create ref_data_singler
colnames(ref_ndata)<- ref_celltype$cell_type
ref_celltype<- ref_celltype$cell_type
ref_maincelltype<- ref_celltype

ref_ndata_singler <- list(name = ref_name, data = as.matrix(ref_ndata),
                         types = ref_celltype, main_types = ref_maincelltype)

ref_ndata_singler$de.genes <- CreateVariableGeneSet(as.matrix(ref_ndata),ref_celltype,200)
ref_ndata_singler$de.genes.main <- CreateVariableGeneSet(as.matrix(ref_ndata),ref_maincelltype,300)

# Annotating
singler_test_data<- CreateSinglerObject(counts = as.matrix(ndata),
                                        project.name = 'Test',
                                        species = species,
                                        ref.list = list(ref_ndata_singler))

singler_test_data <- singler_test_data$singler
singler_test_data <- singler_test_data[[1]]
singler_test_data <- singler_test_data$SingleR.single
celltype$singler <- singler_test_data$labels

