# CellMatch is used to make marker file for each dataset. (https://github.com/ZJUFanLab/scCATCH)
# select a specific species and tissue type for each dataset.
# species ('Human' or 'Mouse'); tissue (e.g.,'Blood','Brain',etc.)

CellMatch<- CellMatch[CellMatch$speciesType == species & CellMatch$tissueType %in% tissue & CellMatch$cancerType == 'Normal',]
cell_types<- unique(CellMatch$cellName)
genename<- unique(CellMatch$geneSymbol)

marker_file<- as.data.frame(matrix(0,nrow = length(genename),ncol = length(cell_types)),stringsAsFactors = F)
rownames(marker_file)<- genename
colnames(marker_file)<- cell_types
for (j in 1:length(cell_types)) {
  d1<- CellMatch[CellMatch$cellName == cell_types[j],]$geneSymbol
  d1<- unique(d1)
  marker_file[d1,j]<- 1
}
