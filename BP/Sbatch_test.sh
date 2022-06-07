#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G  #if use large memory, need lining for some time
#SBATCH --time=0-02:00:00     # 7 days
#SBATCH --mail-user=xliu261@ucr.edu
#SBATCH --mail-type=ALL   #
#SBATCH --job-name="test"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu, short
#SBATCH --output=Sbatch_test.out
#SBATCH --array=1-50  #50 jobs is the best, more jobs won't provide some speed-up



module load itpp
module load anaconda
#echo "${SLURM_ARRAY_TASK_ID}"
python timing_test.py
