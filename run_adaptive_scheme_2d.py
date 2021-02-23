#!/usr/bin/env python3

# to set sgpp paths manually:
# import sys
# sys.path.append('/hppfs/work/pn34mi/di39qun2/SGpp/lib')
# sys.path.append('/hppfs/work/pn34mi/di39qun2/SGpp/lib/pysgpp')
#cf. https://stackoverflow.com/questions/6543847/setting-ld-library-path-from-inside-python
#import os
#os.environ['LD_LIBRARY_PATH']='/hppfs/work/pn34mi/di39qun2/SGpp/lib'
#try:
#    os.execv(sys.argv[0], sys.argv)
#except Exception as exc:
#    print ('Failed re-exec:', exc)
#    sys.exit(1)

import dimadapt2d

# restricted w to 64, because `Gauss-Laguerre integration is currently not possible for nw0>124.`

print("Running with adaptation only to target resolution")
adapt = dimadapt2d.DimensionalAdaptation([5,5,5,5,3], [10,10,10,10,6])

adapt.run_dimadapt_algorithm(100)

