#!/bin/sh

#SBATCH -J fetach_chipseq
#SBATCH --time=5-12:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=12gb

module load bedtools
dir=/mnt/pan/Data14/ASB
#bedfile=${dir}/metafiles/submeta/GM12878_CHiP-seq_default.list${SLURM_ARRAY_TASK_ID}
bedfile=${dir}/metafiles/submeta/asb_missing.list${SLURM_ARRAY_TASK_ID}

bamfile=${dir}/metafiles/GM12878_CHiP-seq_bam_bam.list
bamdir=${dir}/bamfiles/
beddir=${dir}/bedfiles/

cat ${bedfile} | while read line
do
  tf=`echo $line | awk '{print $2}'`
  exp=`echo $line | awk '{print $1}'` 
  bed=`echo $line | awk '{print $NF}'`
  ind=`echo $line | awk '{print $(NF-1)}'`

  bedfile_url=`echo $line | awk '{print $(NF-2)}'`
  wget ${bedfile_url} -O ${beddir}${bed}.gz
  gunzip ${beddir}${bed}.gz

  cat ${bamfile} | grep "$exp $tf" | while read line2
  do
      bam=`echo $line2 | awk '{print $NF}'`
      bamfile_url=`echo $line2 | awk '{print $(NF-2)}'`
      wget ${bamfile_url} -O ${bamdir}${bam}
      samtools index ${bamdir}${bam}
  done

  echo "group_name target replicate_number bam_name bed_name" > ${dir}/samplesheet/${ind}_${exp}_${tf}.csv
  cat ${bamfile} | grep "$exp $tf" | \
	awk -v var="$ind" -v var2="$bamdir" -v var3="$beddir" -v var4="$bed" \
		'{print var,$2,$3,var2$NF,var3""var4}' >> ${dir}/samplesheet/${ind}_${exp}_${tf}.csv
  sed -i 's/\s\+/\t/g' ${dir}/samplesheet/${ind}_${exp}_${tf}.csv
  Rscript lib/allelic_chip-seq.r asb_output/${ind}_${exp}_${tf}_output ${dir}/samplesheet/${ind}_${exp}_${tf}.csv ${dir}/${ind}_hetSNP.txt ${ind}
  rm ${bamdir}*${tf}*.bam 
done

  
