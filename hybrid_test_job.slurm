#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 8
#SBATCH -t 00:10:00
#SBATCH -J hybrid_scaling_test
#SBATCH -e %j.%N.err
#SBATCH -o %j.%N.out
#SBATCH -A m175
#SBATCH --mail-type=ALL,TIME_LIMIT_50
#SBATCH --mail-user=charles@ices.utexas.edu
#SBATCH -V

export MPICH_MAX_THREAD_SAFETY=multiple
export OMP_NUM_THREADS=2
srun -n 96 -N 8  -c 2 ./main > 096-12-02.out
export OMP_NUM_THREADS=4
srun -n 48 -N 8  -c 4 ./main > 048-06-04.out
export OMP_NUM_THREADS=12
srun -n 16 -N 8  -c 12 ./main > 016-02-12.out
export OMP_NUM_THREADS=12
srun -n 8 -N 8  -c 12 ./main > 008-01-12.out
export OMP_NUM_THREADS=24
srun -n 8 -N 8  -c 24 ./main > 008-01-24.out
