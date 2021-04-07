#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# to find the dependencies, run this before
# export PYTHONPATH=$ADAPTATION_GENE3D_PATH/gene-diag/GENE_gui_python:$PYTHONPATH

from utils.loader import Loader
from diag_fluxesgene3d_combi import DiagFluxesgene3dC
from diag_flux_spectragene3d_combi import DiagFluxSpectragene3dC
#from diagnostics.diag_profiles_combi import DiagProfilesC
from diag_profiles_combi import DiagProfilesC
from pathlib import Path
from data.data import Data
from utils.run import Simulation
from commandline import CommandLine

import sys
from os import environ
import numpy as np

#def do_additional_diagnostics(folder, out_folder="./", runextensions=['.dat'], t_start=environ.get("ADAPTATION_PARAM_CROP_TIME"), t_end=float("inf")):
cl=CommandLine(sys.argv[1:])

sim=cl.sim
    
selected_diags = []
selected_diags.append(DiagFluxSpectragene3dC(sim.data.avail_vars,sim.run[0].specnames))
selected_diags.append(DiagFluxesgene3dC(sim.data.avail_vars,sim.run[0].specnames))
selected_diags.append(DiagProfilesC(sim.data.avail_vars,sim.run[0].specnames))

loader = Loader(selected_diags, sim.run[0], sim.data)
loader.set_interval(sim.data, sim.starttime, sim.endtime, sim.stepping)    

##executing
for it, time in enumerate(loader.processable_times):
    print(" time {}".format(time))
    for i_d, diag in enumerate(selected_diags):
        diag.execute(sim.data, loader.processable[it], loader.processable[-1])
    
for i_d, diag in enumerate(selected_diags):
    diag.plot(loader.processable_times, None, cl.sim.out_folder)


## if called directly, run tests
#if __name__ == "__main__":
#    #these are input
#    folder = Path('/hppfs/scratch/02/di39qun2/gene3d-flw-simulations/prob_5_5_7_6_3/out')
#    out_folder = Path('./')
#    runextensions = ['.dat']
#    t_start = 150.
#    t_end = float("inf")
#    do_additional_diagnostics(folder, out_folder, runextensions, t_start, t_end)

