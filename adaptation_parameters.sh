#!/bin/bash

# export any system-wide variables necessary (Python modules, GENE3D paths...)
#TODO


# export adaptation-specific variables
export ADAPTATION_GENE3D_PATH="/hppfs/work/pn34mi/di68xux2/myGene3d/"
export ADAPTATION_PROB_PATH="/hppfs/scratch/02/di39qun2/gene3d-flw-simulations/"

# the parameter file template contains placeholders:
# `$n_procs_x`, `$n_procs_y`, `$n_procs_z`, `$n_procs_v`, `$n_procs_w`, and `$n_procs_s` for the parallelizations
# `$nx0`, `$ny0`, `$nz0`, `$nv0`, and `$nw0` for the discretizations
# `$timelim`, `$simtimelim` for the time limits -- cf `parameters_flw_template` for reference
export ADAPTATION_PARAMETER_FILE_TEMPLATE="/hppfs/work/pn34mi/di39qun2/ParameterFiles/gene3d-qoi/parameters_flw_template"
# the submit script template contains placeholders:
# `$probname`,`$nnodes`,`$nprocesses`,`$p_level`,`$nprocspernode` for the run parameters
# (the parallelization may need to be optimized in `sim_launcher.py`)
# cf `sbatch_template.sh` for reference
export ADAPTATION_SUBMIT_SCRIPT_TEMPLATE="/hppfs/work/pn34mi/di39qun2/ParameterFiles/gene3d-qoi/sbatch_template.sh"
export ADAPTATION_SUBMIT_COMMAND="sbatch"

export ADAPTATION_MINIMUM_LEVEL="[5,5,5,5,3]"