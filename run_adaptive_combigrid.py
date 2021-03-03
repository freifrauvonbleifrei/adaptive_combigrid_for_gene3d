#!/usr/bin/env python3

from os import environ
import dimadapt2d

# restricted w to 64, because `Gauss-Laguerre integration is currently not possible for nw0>124.`

print("Running with adaptation only to target resolution")
adapt = dimadapt2d.DimensionalAdaptation(environ.get('ADAPTATION_MINIMUM_LEVEL'), [10,10,10,10,6], output3d=True)

adapt.run_dimadapt_algorithm(100)
