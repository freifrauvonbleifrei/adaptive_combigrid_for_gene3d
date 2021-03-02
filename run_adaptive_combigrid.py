#!/usr/bin/env python3

import dimadapt2d

# restricted w to 64, because `Gauss-Laguerre integration is currently not possible for nw0>124.`

print("Running with adaptation only to target resolution")
adapt = dimadapt2d.DimensionalAdaptation([5,5,5,5,3], [10,10,10,10,6], output3d=True)

adapt.run_dimadapt_algorithm(100)

