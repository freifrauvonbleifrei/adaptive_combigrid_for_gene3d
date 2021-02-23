#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J $probname
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
##SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=theresa.pollinger@ipvs.uni-stuttgart.de
# Wall clock limit:
#SBATCH --time=20:00:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn34mi

##SBATCH --partition=fat
#SBATCH --partition=micro
##SBATCH --partition=test
#constraints are optional
#--constraint="scratch&work"
#Number of nodes and MPI tasks per node:
#SBATCH --nodes=$nnodes
#SBATCH --ntasks=$nprocesses
##SBATCH --ntasks-per-node=$nprocspernode
#SBATCH --cpus-per-task=1

umask 002
 
#Important
module load slurm_setup
module switch mpi.intel mpi.intel/2019.7
module load hdf5
 
#Run the program:
mpiexec -n $SLURM_NTASKS ./gene3d -no_signal_handler
