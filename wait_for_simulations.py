#!/usr/bin/env python3

import sys
import sim_launcher
import time
from Qes_data import Qes_data
from write_flux_adaptive_combigrid import get_combiScheme

combiSchemeMode = sys.argv[1]
newEndTime = sys.argv[2]

combiScheme = get_combiScheme("", combiSchemeMode, dropzeros=True)

d = 0
while str('l_'+str(d+1)) in combiScheme.columns:
    d += 1

while len(combiScheme) > 0:
    for i, row in combiScheme.iterrows():
        level = [0]*d
        lkeys = [str('l_'+str(dim)) for dim in range(d)]
        for dim in range(d):
            level[dim] = combiScheme[lkeys[dim]]
        # check if current simulation is finished
        if sim_launcher.check_finished(level_vector=level):
            # if yes, check if simulation time is reached
            if Qes_data(level_vector=level).how_long_run() > 0.99 * newStartTime:
                # if yes, remove from list
                combiScheme.drop(index=i, inplace=True)
            else:
                # if no, prolong the simulation again
                sim_launcher.restart_sim(
                    level_vector=level, newSimTime=newEndTime)
    # sleep for a minute, so we don't burn energy absolutely unnecessarily
    time.sleep(60)
