#!/bin/sh

wget https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv
cat tabix_ftp_paths.tsv  | grep -e "LCL" | grep -e "naive" | grep -e "ge" | awk -F"\t" '{print $NF}' | while read file
do 
  wget $file 
done

ls *.tsv.gz |  while read file
do
zcat $file | awk '$9<10e-6 || (NR==1){print}' | grep -e "SNP"  -e "gene_id" > ${file%all.tsv.gz}thres.tsv
done


