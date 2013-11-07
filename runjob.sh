#!/bin/bash
# The name of the job.
#PBS -N test
# Asking for 12 cores.
#PBS -l nodes=12
# Mail user if job aborts (a) or ends (e) (currently commented out).
# PBS -m ae
# Set the job to run 24 hours 10 minutes and 5 seconds.
#PBS -l walltime=00:05:0
# Set output file name.
#PBS -o test.txt
#PBS -V
# module load valgrind/3.8.1-openmpi-gcc
module load gcc/4.8.1
module load openmpi/1.7.3-gcc
module unload gcc/4.7.2
cd $PBS_O_WORKDIR
# LD_LIBRARY_PATH=$PBS_O_WORKDIR/spatialindex/lib:$LD_LIBRARY_PATH mpiexec ./particle -l -10 -r 10 -d 10 -n 174
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PBS_O_WORKDIR/spatialindex/lib mpiexec ./particle -n 10 -s ROSENBR
