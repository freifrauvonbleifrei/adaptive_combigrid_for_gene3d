#!/usr/bin/env python3

from Qes_data import Qes_data
import sim_launcher
from write_flux_adaptive_combigrid import get_combiScheme
import sys


def prolong(combiSchemeMode, newStartTime, newEndTime, checkTime=True):
    combiScheme = get_combiScheme("", combiSchemeMode, dropzeros=True)

    d = 0
    while str('l_'+str(d)) in combiScheme:
        d += 1

    level = [0]*d
    lkeys = [str('l_'+str(dim)) for dim in range(d)]

    for i, row in combiScheme.iterrows():
        for dim in range(d):
            level[dim] = row[lkeys[dim]]
        # check if current simulation is finished
        print("trying to prolong " + str(level))
        assert(sim_launcher.check_finished(level_vector=level))
        if checkTime:
            # and has run long enough
            assert(Qes_data(level_vector=level).how_long_run() > 0.99 * newStartTime)
            # but not too long
            assert(Qes_data(level_vector=level).how_long_run() < 1.01 * newStartTime)
        sim_launcher.restart_sim(level_vector=level, newSimTime=newEndTime)

if __name__ == "__main__":
    combiSchemeMode = sys.argv[1]
    newStartTime = sys.argv[2]
    newEndTime = sys.argv[3]
    prolong(combiSchemeMode, newStartTime, newEndTime)
