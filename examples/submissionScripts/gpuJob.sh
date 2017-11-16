#!/bin/bash

# set the number of nodes and processes per node. We are running one process on a single node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# set max wallclock time
#SBATCH --time=00:30:00

# set name of job
#SBATCH --job-name QUEST

# set queue
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
###SBATCH --partition=devel-gpu

NUM_QUBITS=27
EXE=demo
export OMP_NUM_THREADS=16

module purge
module load gpu/cuda

## uncomment for reporting about GPU used
##/system/software/arcus-b/gpu/cuda/7.0.28/samples/1_Utilities/deviceQuery/deviceQuery

make clean
make

./$EXE $NUM_QUBITS

## uncomment to run a command line profiler
##nvprof ./$EXE $NUM_QUBITS
