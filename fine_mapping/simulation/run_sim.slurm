#!/bin/sh

#SBATCH -J simulation
#SBATCH --time=5-12:00:00
#SBATCH -N 1 -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 6g

maf=0.2
for m in $(seq 100 100 1000)
do
  for n in $(seq 8 2 30)
  do
    python lib/setup1.py ${maf} ${n} ${m} kp1.list conf1_kp1_cell${m}_sub${n}_${SLURM_ARRAY_TASK_ID}
    sed -i "s/$/,${SLURM_ARRAY_TASK_ID}/g" conf1_kp1_cell${m}_sub${n}_${SLURM_ARRAY_TASK_ID}
#    python lib/setup2.py ${maf} ${n} ${m} kp1.list conf2_kp1_cell${m}_sub${n}_${SLURM_ARRAY_TASK_ID}
#    sed -i "s/$/,${SLURM_ARRAY_TASK_ID}/g" conf2_kp1_cell${m}_sub${n}_${SLURM_ARRAY_TASK_ID}
#    python lib/setup3.py ${maf} ${n} ${m} kp1.list conf3_kp1_cell${m}_sub${n}_${SLURM_ARRAY_TASK_ID}
#    sed -i "s/$/,${SLURM_ARRAY_TASK_ID}/g" conf3_kp1_cell${m}_sub${n}_${SLURM_ARRAY_TASK_ID}
  done
done 

#cat conf1_kp1_cell* > conf1_kp1.res
#cat conf2_kp1_cell* > conf2_kp1.res
#cat conf3_kp1_cell* > conf3_kp1.res

#python scatter.py conf1_kp1.res 0.05 conf1_kp1
#python scatter.py conf1_kp1.res 1e-6 conf1_kp1

#python scatter.py conf2_kp1.res 0.05 conf2_kp1
#python scatter.py conf2_kp1.res 1e-6 conf2_kp1










 
