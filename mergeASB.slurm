#!/bin/sh

#SBATCH -J mergeASB
#SBATCH --time=1-12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=6gb

module load bedtools
dir=/mnt/pan/Data14/ASB
#bedfile=${dir}/metafiles/submeta/GM12878_CHiP-seq_default.list${SLURM_ARRAY_TASK_ID}
bedfile=${dir}/metafiles/submeta/asb_missing.list${SLURM_ARRAY_TASK_ID}
bamdir=$pan/ASB/bamfiles/
beddir=$pan/ASB/bedfiles/
 
cat ${bedfile} | while read line
do
  tf=`echo $line | awk '{print $2}'`
  exp=`echo $line | awk '{print $1}'`
  bed=`echo $line | awk '{print $NF}'`
  ind=`echo $line | awk '{print $(NF-1)}'`

  asb_file=asb_output/${ind}_${exp}_${tf}_output.res
  bed_file=bedfiles/${ind}_${exp}-${tf}.bed
  vcf=/mnt/pan/Data14/1kgp/NA${ind#GM}.GRCh38.phased.vcf.het.gz
  python lib/merge_asb.py mergePeak ${asb_file} ${ind}_hetSNP.bed ${bed_file} ${ind} ${asb_file%.res}
  sed -i "s/${ind}\.//g" ${asb_file%.res}.merge
  python lib/merge_asb.py filterASB ${asb_file%.res}.merge ${vcf} ${dir}/asb_output/${ind}_${exp}_${tf}
done


