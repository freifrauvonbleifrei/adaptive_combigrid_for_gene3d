#!/usr/bin/env python3

from os import environ
import dimadapt2d

# restricted w to 64, because `Gauss-Laguerre integration is currently not possible for nw0>124.`

print("Running with adaptation only to target resolution")
# cf. https://stackoverflow.com/questions/48066106/converting-vector-looking-string-into-vector-values-python
lmin=[int(i.strip()) for i in environ.get('ADAPTATION_MINIMUM_LEVEL')[1:-1].split(",")]
lmax=[int(i.strip()) for i in environ.get('ADAPTATION_MAXIMUM_LEVEL')[1:-1].split(",")]
adapt = dimadapt2d.DimensionalAdaptation(lmin, lmax, output3d=True)

adapt.run_dimadapt_algorithm(int(environ.get('ADAPTATION_PARAM_MAXIMUM_NUM_GRIDS')))
