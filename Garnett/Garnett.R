library(garnett)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# ndata_cds is obtained from [cds_object_making.R]
# marker_file.txt is obtained from [marker_file_making.R]
# data_db (org.Hs.eg.db for human dataset; org.Mm.eg.db for mouse datasets)

garnett_classifier <- train_cell_classifier(cds = ndata_cds,
                                            marker_file = marker_file.txt,
                                            db = data_db,
                                            cds_gene_id_type = "SYMBOL",
                                            num_unknown = 50,
                                            marker_file_gene_id_type = "SYMBOL")

ndata_cds <- classify_cells(ndata_cds, garnett_classifier,
                            db = data_db,
                            cluster_extend = TRUE,
                            cds_gene_id_type = "SYMBOL")
