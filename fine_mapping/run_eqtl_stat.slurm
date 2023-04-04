#!/bin/sh

#SBATCH -J eqtl_corr
#SBATCH --time=1-12:00:00
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem=30gb

sample=$1
dir=/mnt/pan/Data14
python eqtl_stat.py all_LCL.thres.tsv ${dir}/LCL_kinetics/${sample}.GRCh38.phased.het.vcf.gz \
		${dir}/LCL_kinetics/${sample}/${sample}_pb_qc.est ${sample}_all_LCL.thres 

python kp_stat.py ${dir}/LCL_kinetics/NA12878/NA12878_pb_qc.est,${dir}/LCL_kinetics/NA18502/NA18502_pb_qc.est NA12878,NA18502 LCL 
