#!/bin/bash
#SBATCH --output=version1_%j.out
#SBATCH --job-name=remove_g_v1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

start_time=$(date +%s)

echo "Job started at: $(date)"

python3 remove_g_updated.py $1 $2

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "Job finished at: $(date)"
echo "Total runtime: ${runtime} seconds"
