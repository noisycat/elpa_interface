#!/bin/bash
#PBS -q debug
#PBS -l mppwidth=48
#PBS -l walltime=00:30:00
#PBS -N Si525
#PBS -e $PBS_JOBID.err
#PBS -o $PBS_JOBID.out
#PBS -A m175
#PBS -V

cd $PBS_O_WORKDIR
export MPICH_MAX_THREAD_SAFETY multiple
export OMP_NUM_THREADS=24
aprun -n 4 -N 2 -ss -d 12 valgrind --leak-check=yes ./main 8000
