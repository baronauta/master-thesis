#!/bin/bash
#SBATCH --job-name=test
#SBATCH --error=test_%j.err
#SBATCH --output=test_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=24:00:00

cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}

MKL_NUM_THREADS=2 julia myscript

echo "calculation finished"
exit 0