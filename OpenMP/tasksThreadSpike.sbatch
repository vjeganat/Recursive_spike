#!/bin/bash -l
#SBATCH -N 1        # you'll never need more than 1 node for the openmp code
#SBATCH -n 20
#SBATCH -t 00:05:00   # adjust the amount of time as necessary
#SBATCH -p gpu
#SBATCH -A UH-COSC6374-18F-Para
#SBATCH -J openmp
#SBATCH -o TasksOmp.o%j


export OMP_PLACES=cores
export OMP_PROC_BIND=close

./spike 65536 1
./spike 65536 2
./spike 65536 4
./spike 65536 8
./spike 65536 16
./spike 65536 32
./spike 65536 64
./spike 65536 128
./spike 65536 256
./spike 65536 512



