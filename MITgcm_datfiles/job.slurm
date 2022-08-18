#!/bin/bash
#
# Output directory
#SBATCH -o /data/directory/.out
#SBACTH -e data/directory/.err
#SBATCH -D data/directory/
#
#SBATCH --mem=16gb # Maximum memory 
#SBATCH --time=24:00:00 # Maximum CPU time
#SBATCH -n 100 # number of tasks
#SBATCH -N 2  # number of nodes
##SBATCH --mem-per-cpu=2gb
#
# Job name
#SBATCH -J job-name
#
# Send mail when the job beings, ends, fails or is requeued
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=email@email.co.uk
#
# Pick the partition (equivalent to SGE queue) to run the job in.
#SBATCH --partition=medium
#SBATCH --account=medium

# Load the modules
source /etc/profile.d/modules.sh
module load hpc/mvapich2/gcc/2.2 hpc/netcdf/gcc/4.4.1.1 hpc/openmpi/gcc/3.0.0 hpc/gcc/7.2.0

# unlimit the stack to made job run.
ulimit -s hard

# Run the job on 5 cores.
srun -n100 --mpi=pmi2 ./mitgcmuv
