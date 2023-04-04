#!/bin/bash

sample=$1
dir1=/mnt/pan/Data14/LCL_kinetics
dir2=/mnt/pan/Data14/rosmap_microglia
i=1
cut -d, -f1 ${dir1}/${sample}/outs/aggregation.csv | sed 1d | while read id
do
  if [[ $id =~ "Microglia" ]]
  then
    sed "s/-1/-${i}/" ${dir2}/sc_${id}/${id}_qc.csv >> ${dir1}/${sample}/${sample}_qc.csv
    sed "s/-1/-${i}/" ${dir2}/sc_${id}/${id}_qc.csv | cut -d, -f1 | \
	awk '{print $0",Microglia"}' >> ${dir1}/${sample}/${sample}_qc.ct
  else
    sed "s/-1/-${i}/" ${dir1}/${id}/${id}_qc.csv >> ${dir1}/${sample}/${sample}_qc.csv
    sed "s/-1/-${i}/" ${dir1}/${id}/${id}_qc.csv | cut -d, -f1 | \
	awk '{print $0",LCL"}' >> ${dir1}/${sample}/${sample}_qc.ct
  fi
  i=`expr $i + 1` 
done

sed -i '1!{/cb,/d}' ${dir1}/${sample}/${sample}_qc.csv
sed -i '1!{/cb,/d}' ${dir1}/${sample}/${sample}_qc.ct

