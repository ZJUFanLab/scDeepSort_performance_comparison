library(garnett)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# CellMatch is used to make marker file for each dataset. (https://github.com/ZJUFanLab/scCATCH)
# select a specific species and tissue type for each dataset.
# species ('Human' or 'Mouse'); tissue (e.g.,'Blood','Brain',etc.)
CellMatch<- CellMatch[CellMatch$speciesType == species & CellMatch$tissueType %in% tissue & CellMatch$cancerType == 'Normal',]
cell_types<- unique(CellMatch$cellName)

marker_file.txt <- NULL
for (j in 1:length(cell_types)) {
  marker_file.txt<- c(marker_file.txt,paste('>',cell_types[j],sep = ''))
  
  d1<- CellMatch[CellMatch$cellName == cell_types[j],]$geneSymbol
  d1<- unique(d1)
  res_markers<- NULL
  if (length(d1) == 1) {
    res_markers<- d1
  }
  if (length(d1) > 1) {
    res_markers<- d1[1]
    for (k in 2:length(d1)) {
      res_markers<- paste(res_markers,d1[k],sep = ', ')
    }
  }
  res_markers<- paste('expressed: ',res_markers,sep = '')
  marker_file.txt<- c(marker_file.txt,res_markers)
  
  d2<- CellMatch[CellMatch$cellName == cell_types[j],]$PMID
  d2<- unique(d2)
  res_reference<- NULL
  if (length(d2) == 1) {
    res_reference<- d2
  }
  if (length(d2) > 1) {
    res_reference<- d2[1]
    for (k in 2:length(d2)) {
      res_reference<- paste(res_reference,d2[k],sep = ', ')
    }
  }
  res_reference<- paste('references: ',res_reference,sep = '')
  marker_file.txt<- c(marker_file.txt,res_reference,'')
}

# ndata_cds is obtained from [cds_object_making.R]
# check and reconstruct marker file for each dataset.
# data_db (org.Hs.eg.db for human dataset; org.Mm.eg.db for mouse datasets)
marker_check <- check_markers(ndata_cds, marker_file.txt,
                              db = data_db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

# exluding markers of Low nomination, High ambiguity, Not in db, Not in CDS.
marker_check$in_cds<- as.character(marker_check$in_cds)
marker_check<- marker_check[marker_check$summary == 'Low nomination?' | marker_check$summary == 'High ambiguity?' | marker_check$summary == 'Not in db' | marker_check$summary == 'Not in CDS',]
CellMatch<- CellMatch[!CellMatch$geneSymbol %in% marker_check$marker_gene,]

# re-making marker file

marker_file.txt <- NULL
for (j in 1:length(cell_types)) {
  marker_file.txt<- c(marker_file.txt,paste('>',cell_types[j],sep = ''))
  
  d1<- CellMatch[CellMatch$cellName == cell_types[j],]$geneSymbol
  d1<- unique(d1)
  res_markers<- NULL
  if (length(d1) == 1) {
    res_markers<- d1
  }
  if (length(d1) > 1) {
    res_markers<- d1[1]
    for (k in 2:length(d1)) {
      res_markers<- paste(res_markers,d1[k],sep = ', ')
    }
  }
  res_markers<- paste('expressed: ',res_markers,sep = '')
  marker_file.txt<- c(marker_file.txt,res_markers)
  
  d2<- CellMatch[CellMatch$cellName == cell_types[j],]$PMID
  d2<- unique(d2)
  res_reference<- NULL
  if (length(d2) == 1) {
    res_reference<- d2
  }
  if (length(d2) > 1) {
    res_reference<- d2[1]
    for (k in 2:length(d2)) {
      res_reference<- paste(res_reference,d2[k],sep = ', ')
    }
  }
  res_reference<- paste('references: ',res_reference,sep = '')
  marker_file.txt<- c(marker_file.txt,res_reference,'')
}
