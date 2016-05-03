#!/bin/bash


split_libraries.py \
	-m 040816JP27F-mapping.txt \
	-f data_raw/fasta-qual-mapping-files/040816JP27F-full.fasta \
	-q data_raw/fasta-qual-mapping-files/040816JP27F-full.qual -o split_library_output -b 8
