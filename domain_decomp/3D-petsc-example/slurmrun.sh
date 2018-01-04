#!/bin/bash
export jobName=test
#Account and Email Information
#SBATCH --mail-type=end
#SBATCH --mail-user=scott@aiton.net

# Specify parition (queue)
#SBATCH --partition=batch

# Join output and errors into output
#SBATCH -o cout.o%j
#SBATCH -e cerr.e%j

# Specify job not to be rerunable
#SBATCH --no-requeue

# Specify walltime
###SBATCH --time=48:00:00

# Specify the procs per node:

#Exclusively check out a node
#SBATCH --exclusive

#Module Load Section
module load openmpi/2.0.1
module load cuda75/7.5
module load gcc/4.8.1 
module load extra
module load my-petsc

cd $SLURM_SUBMIT_DIR

export CUDA_LAUNCH_BLOCKING=1

srun --mpi=pmi2 $@
sleep 10s
