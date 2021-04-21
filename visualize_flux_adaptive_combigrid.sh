#!/bin/bash

. ./spack-sgpp/share/spack/setup-env.sh
spack load sgpp

. ./adaptation_parameters.sh

python3 visualize_flux_adaptive_combigrid.py $1
