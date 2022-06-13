#!/bin/bash -l
#SBATCH -N 1        # you'll never need more than 1 node for the openmp code
#SBATCH -n 20
#SBATCH -t 00:05:00   # adjust the amount of time as necessary
#SBATCH -p gpu
#SBATCH -A UH-COSC6374-18F-Para
#SBATCH -J openmp
#SBATCH -o dynamicOmp.o%j


export OMP_PLACES=cores
export OMP_PROC_BIND=close

./dspike 512
./dspike 4096
./dspike 65536
./dspike 16777216



