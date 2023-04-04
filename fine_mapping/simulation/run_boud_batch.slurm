#!/bin/sh

#SBATCH -J simulation
#SBATCH --time=5-12:00:00
#SBATCH -N 1 -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 6g

cell=$1
sub=$2

for logfc_bs in $(seq -2 0.5 2)
do
  for logfc_bf in $(seq -2 0.5 2)
  do
   python lib/boundary_simulation.py kp1.list ${logfc_bs} ${logfc_bf} ${cell} ${sub}\
	bond_cell${cell}_sub${sub}_kp1_bs${logfc_bs}_bf${logfc_bf}_${SLURM_ARRAY_TASK_ID}
   sed -i "s/$/,${SLURM_ARRAY_TASK_ID}/g" \
	bond_cell${cell}_sub${sub}_kp1_bs${logfc_bs}_bf${logfc_bf}_${SLURM_ARRAY_TASK_ID} 
 done
done 

 
