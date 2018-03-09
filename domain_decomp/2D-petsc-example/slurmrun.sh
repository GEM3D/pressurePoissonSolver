#!/bin/bash
export jobName=test
#Account and Email Information

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

cd $SLURM_SUBMIT_DIR

srun --mpi=pmi2 $@
sleep 10s
