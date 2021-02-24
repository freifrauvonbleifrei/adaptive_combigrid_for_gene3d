#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J adaptive_update 
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
##SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=theresa.pollinger@ipvs.uni-stuttgart.de
# Wall clock limit:
#SBATCH --time=00:05:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn34mi

#SBATCH --partition=test
#constraints are optional
#--constraint="scratch&work"
#Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

umask 002
 
#Important
module load slurm_setup
module load python/3.6_intel

. ./spack-sgpp/share/spack/setup-env.sh
spack load sgpp

. ./adaptation_parameters.sh

python3 run_adaptive_scheme_2d.py
