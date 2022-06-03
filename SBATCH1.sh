#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=256G
#SBATCH --time=0-02:00:00     # 7 days
#SBATCH --mail-user=xliu261@ucr.edu
#SBATCH --mail-type=ALL   #
#SBATCH --job-name="pro0110_12time"
#SBATCH -p short # This is the default partition, you can use any of the following; intel, batch, highmem, gpu
#SBATCH --output=pro0110_12time.out


# Print current date


# Load samtools
module load itpp

g++  -std=c++11 -fopenmp `itpp-config --cflags` -o my_pro product.cpp `itpp-config --libs`

./my_pro




# Concatenate BAMs
#samtools cat -h header.sam -o out.bam in1.bam in2.bam

# Print name of node
