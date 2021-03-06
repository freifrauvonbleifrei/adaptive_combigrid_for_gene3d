#!/bin/bash

# export any system-wide variables necessary (Python modules, GENE3D paths...)
#TODO


# export adaptation-specific variables
export ADAPTATION_GENE3D_PATH="/hppfs/work/pn34mi/di68xux2/myGene3d/"
export ADAPTATION_PROB_PATH="/hppfs/scratch/02/di39qun2/gene3d-flw-simulations/"

# specify the number of species (currently only 0 for adiabatic electrons and 1 for kinetic electrons allowed)
export ADAPTATION_NUMBER_OF_SPECIES="2"
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
export ADAPTATION_MAXIMUM_LEVEL="[8,8,9,9,6]"

export ADAPTATION_PARAM_Y_EQUALS_X=true
export ADAPTATION_PARAM_MAXIMUM_NUM_GRIDS="20"
export ADAPTATION_PARAM_CROP_TIME="150.0"
export ADAPTATION_PARAM_MIN_RUNTIME="500.0"
export ADAPTATION_PARAM_TERMINATION_CRITERION_SEM_FACTOR="1."
export ADAPTATION_PARAM_SEM_MINIMUM_NUMBER_WINDOWS="3"
# relative adaptation: the delta for adaptation will be divided by the component grid's value 
# thus favoring partial solutions with lower values
export ADAPTATION_PARAM_RELATIVE_ADAPTATION=true
# relative combination: the partial solutions for combination will be normalized to the total qes value
# thus avoiding to add high relative noise in the combined solution
export ADAPTATION_POSTPROCESSING_RELATIVE_COMBINATION=true
# rolling average: smooth the profile with its neighbors before (rescaling and) combining
# =1 corresponds to no averaging
export ADAPTATION_POSTPROCESSING_ROLLING_AVG_NUM_POINTS=1

export ADAPTATION_RESULTS_CSV="qes_results_species.csv"
