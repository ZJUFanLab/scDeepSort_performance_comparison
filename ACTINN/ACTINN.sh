#!/bin/sh

python actinn_format.py -i ref_ndata_train.csv.gz -o ref_ndata_train -f csv
python actinn_format.py -i ndata.csv.gz -o ndata -f csv
python actinn_predict.py -trs ref_ndata_train.h5 -trl ref_celltype.txt.gz -ts ndata.h5 -lr 0.0001 -ne 50 -ms 128 -pc True -op False