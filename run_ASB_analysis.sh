#!/bin/sh

dir=/mnt/pan/Data14/ASB
ind=GM12878
bedfile=${dir}/metafiles/${ind}_CHiP-seq_default.list
bamdir=$pan/ASB/bamfiles/
beddir=$pan/ASB/bedfiles/

header=`cat ${dir}/asb_output/${ind}*_isASB.final | head -1` 
echo $header",TF" > ${dir}/${ind}_regulator_isASB.final
cat ${bedfile} | while read line
do
  tf=`echo $line | awk '{print $2}'`
  exp=`echo $line | awk '{print $1}'`
  bed=`echo $line | awk '{print $NF}'`
  ind=`echo $line | awk '{print $(NF-1)}'`
  cat ${dir}/asb_output/${ind}_${exp}_${tf}_isASB.final | sed 1d | awk -v var="$tf" '{print $0","var}' >> ${dir}/${ind}_regulator_isASB.final 
done

estimation=/mnt/pan/Data14/LCL_kinetics/NA12878/NA12878_pb_qc.est
python -W ignore lib/analysis_pair.py GM12878_regulator_isASB.final qbic_pred/GM12878_distalEP_hg38.vcf ${estimation} results_v2/GM12878_distalEP_regulator
python -W ignore lib/analysis_pair.py GM12878_regulator_isASB.final qbic_pred/GM12878_TSS2k_hg38.vcf ${estimation} results_v2/GM12878_TSS2k_regulator
python -W ignore lib/lmm.py GM12878_regulator_isASB.final qbic_pred/GM12878_TSS2k-distalEP_hg38.vcf ${estimation} results_v2/GM12878_TSS2k-distalEP_regulator
python lib/barplot_h2.py results_v2/GM12878_TSS2k_regulator_h2.csv results_v2/GM12878_distalEP_regulator_h2.csv

#python -W ignore lib/analysis_ASM.py GM12878_ASM_hg38CG_hg38SNP.bed qbic_pred/GM12878_TSS2k-distalEP_hg38.vcf ${estimation} results/GM12878_TSS2k-distalEP_met

python lib/ASB_eqtl.py GM12878_regulator_isASB.final fine_mapping/NA12878_all_LCL.thres  \
	fine_mapping/NA12878_GC.tsv  GM12878 qbic_pred/GM12878_phased_het_qbic.res qbic_pred/GM18502_dist50k_hg38.vcf
python lib/ASM_eqtl.py GM12878_ASM_hg38CG_hg38SNP.bed fine_mapping/NA12878_all_LCL.thres  fine_mapping/NA12878_GC.tsv  GM12878

python lib/dkp_analysis.py /mnt/pan/Data14/LCL_kinetics/NA12878/NA12878_pb_qc.est,/mnt/pan/Data14/LCL_kinetics/NA18502/NA18502_pb_qc.est \
		NA12878,NA18502 ENSG00000196735 LCL_kinetics_hist
