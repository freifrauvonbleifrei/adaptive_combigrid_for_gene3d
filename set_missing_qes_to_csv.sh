#!/bin/bash

. ./spack-sgpp/share/spack/setup-env.sh
spack load sgpp

. ./adaptation_parameters.sh

python3 set_missing_qes_to_csv.py