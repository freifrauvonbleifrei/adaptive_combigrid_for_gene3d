#!/bin/bash

if [ $# -ne 2 ]; then
    echo "You need two arguments: the new start (1) and end time (2)"
    exit 1
fi

. ./spack-sgpp/share/spack/setup-env.sh
spack load sgpp

. ./adaptation_parameters.sh

python3 prolong_simulations.py oldSet.csv $1 $2

python3 wait_for_simulations.py oldSet.csv $2

python3 write_flux_adaptive_combigrid.py oldSet.csv $1 $2
