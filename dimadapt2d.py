import os
import numpy as np
import pandas as pd
#import seaborn as sns
import math
# import kinds

# this needs `spack load sgpp@adaptive_combigrid_convenience` 
import pysgpp

import Qes_data
import sim_launcher

# utilities for dealing with swig-wrapped c++-vectors
def thingToString(thing):
    s = ""
    for t in thing:
        s += str(t) + "; "
    return s

def thingToList(thing, output=False):
    l = []
    for i in range(len(thing)):
        if output:
            print(thing[i])
        l.append(thing[i])
    return l

def thingToStringList(thing):
    l = []
    for t in thing:
        l.append(str(t))
    return l


class DimensionalAdaptation:
    def __init__(self, lmin, lmax, omega=1., \
                prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'),\
                gene_path=os.environ.get('ADAPTATION_GENE3D_PATH'),\
                output3d=False):
        self.prob_prepath=prob_prepath
        self.gene_path=gene_path
        self.output3d=output3d
        self.lmin=lmin
        self.lmax=lmax
        self.minimum_time_length=float(os.environ.get('ADAPTATION_PARAM_MIN_RUNTIME'))
        #self.minimum_time_length=500
        self.termination_criterion_sem_factor=float(os.environ.get('ADAPTATION_PARAM_TERMINATION_CRITERION_SEM_FACTOR'))
        #self.termination_criterion_sem_factor=1.5
        self.minimum_nwin=int(os.environ.get('ADAPTATION_PARAM_SEM_MINIMUM_NUMBER_WINDOWS'))
        #self.minimum_nwin=3

        dim = len(lmin)
        lmin = pysgpp.LevelVector(lmin)
        levels = pysgpp.LevelVectorVector(1,lmin)
        weightedRelevanceCalculator=pysgpp.WeightedRelevanceCalculator(omega)
        #TODO allow passing relevance generator with own omega -- for now setting omega to 1. in sgpp manually
        assert(omega == 1.)
        self.adaptiveGeneratorElectrons = pysgpp.AdaptiveCombinationGridGenerator(levels)#, weightedRelevanceCalculator)
        self.adaptiveGeneratorIons = pysgpp.AdaptiveCombinationGridGenerator(levels)#, weightedRelevanceCalculator)
        
        # set start results
        qes_data = Qes_data.Qes_data(lmin)
        result = qes_data.get_result()
        sgppActiveLevelVector=pysgpp.LevelVector(lmin)
        self.adaptiveGeneratorElectrons.setQoIInformation(sgppActiveLevelVector, result[1])
        self.adaptiveGeneratorIons.setQoIInformation(sgppActiveLevelVector, result[0])
        
        #same for SEMs
        sem = qes_data.get_time_error()
        self.adaptiveSEMElectrons = \
            pysgpp.AdaptiveCombinationGridGenerator.fromFunctionPointer(levels, \
                                        pysgpp.DoubleVector(1,sem[1]), pysgpp.squaredSummation_cb)
        self.adaptiveSEMIons = pysgpp.AdaptiveCombinationGridGenerator.fromFunctionPointer(levels, \
                                        pysgpp.DoubleVector(1,sem[0]), pysgpp.squaredSummation_cb)
        delta = [self.adaptiveGeneratorElectrons.getDelta(sgppActiveLevelVector),
                self.adaptiveGeneratorIons.getDelta(sgppActiveLevelVector)]
        qes_data.set_delta_to_csv(delta)

        print("starting with: " + str(self.adaptiveGeneratorElectrons.getCurrentResult()) \
                       + "; "  + str(self.adaptiveGeneratorIons.getCurrentResult())\
                      + "] +- [" + str(self.adaptiveSEMElectrons.getCurrentResult()) \
                      + "; "  + str(self.adaptiveSEMIons.getCurrentResult()) + "]"\
        )

        if ((qes_data.how_long_run() < self.minimum_time_length) or qes_data.how_many_nwin_run() < self.minimum_nwin):
            print("please run the initial simulation longer!")
            assert(False)


    def run_dimadapt_algorithm(self, numGrids=10):
        absAdaptedDelta = 100000. # kinds.default_float_kind.MAX
        totalSEM=-1.
        waitingForResults=False
        waitingForNumberOfResults=0
        while (len(self.adaptiveGeneratorElectrons.getOldSet()) < numGrids) and \
                (totalSEM < self.termination_criterion_sem_factor* absAdaptedDelta) and \
                waitingForNumberOfResults < 5 and \
                ~(len(self.adaptiveGeneratorElectrons.getRelevanceOfActiveSet()) == 0):
            waitingForResults=True
            waitingForNumberOfResults=0            

            # iterate the whole active set
            for activeLevelVector in self.adaptiveGeneratorElectrons.getActiveSet():
                #print(activeLevelVector, waitingForNumberOfResults)
                #print(not (self.adaptiveGeneratorElectrons.hasQoIInformation(pysgpp.LevelVector(activeLevelVector))))
                # store QoI information if we have not done it already 
                
                if (not self.adaptiveGeneratorElectrons.hasQoIInformation(pysgpp.LevelVector(activeLevelVector))):
                    print(activeLevelVector, waitingForNumberOfResults)
                    qes_data = Qes_data.Qes_data(activeLevelVector)
                    result = qes_data.get_result()
    #                result = [0.3,0.3]
                    if result:
                        waitingForResults=False
                        sgppActiveLevelVector=pysgpp.LevelVector(activeLevelVector)
                        # store the results in the adaptive generators
                        self.adaptiveGeneratorElectrons.setQoIInformation(sgppActiveLevelVector, result[1])
                        self.adaptiveGeneratorIons.setQoIInformation(sgppActiveLevelVector, result[0])
                        sem=qes_data.get_time_error()
                        self.adaptiveSEMElectrons.setQoIInformation(sgppActiveLevelVector, sem[1])
                        self.adaptiveSEMIons.setQoIInformation(sgppActiveLevelVector, sem[0])
                        print("has values: " + str(result[0]) \
                            + "; "  + str(result[1])\
                            + "] +- [" + str(sem[0]) \
                            + "; "  + str(sem[1]) + "]"\
                        )
                        delta = [self.adaptiveGeneratorElectrons.getDelta(sgppActiveLevelVector),
                                self.adaptiveGeneratorIons.getDelta(sgppActiveLevelVector)]
                        qes_data.set_delta_to_csv(delta)
                        if ((qes_data.how_long_run() < self.minimum_time_length) or qes_data.how_many_nwin_run() < self.minimum_nwin):
                        #prolong simulation
                            if sim_launcher.check_folder_exists(self.prob_prepath, activeLevelVector):
                                if sim_launcher.check_finished(self.prob_prepath, activeLevelVector):
                                    # restart the simulation
                                    print("sim_launcher: restart "+str(activeLevelVector))
                                    sim_launcher.restart_sim(self.prob_prepath, self.gene_path, activeLevelVector)
                        elif (sem[0] > 0.07*result[0]) or (sem[1] > 0.07*result[1]):
                        #prolong simulation
                            if sim_launcher.check_folder_exists(self.prob_prepath, activeLevelVector):
                                if sim_launcher.check_finished(self.prob_prepath, activeLevelVector):
                                    # restart the simulation
                                    print("sim_launcher: restart because of SEM "+str(activeLevelVector))
                                    sim_launcher.restart_sim(self.prob_prepath, self.gene_path, activeLevelVector, qes_data.how_long_run() + self.minimum_time_length)

                    else:
                        # start or prolong simulations
                        if sim_launcher.check_folder_exists(self.prob_prepath, activeLevelVector):
                            if sim_launcher.check_finished(self.prob_prepath, activeLevelVector):
                                # restart the simulation
                                print("sim_launcher: restart "+str(activeLevelVector))
                                sim_launcher.restart_sim(self.prob_prepath, self.gene_path, activeLevelVector)
                            # else, we are just waiting for the simulation to run long enough
                        else:
                            maxres = False
                            for a,l in zip(activeLevelVector, self.lmax):
                                if a > l:
                                    maxres = True
                            if maxres:
                                print("reached maximum resolution at "+str(activeLevelVector))
                            else:
                                # start a simulation
                                print("sim_launcher: start "+str(activeLevelVector))
                                sim_launcher.dispatch_to_run(self.prob_prepath, activeLevelVector, self.gene_path)
                        waitingForNumberOfResults+=1
                else: 
                    waitingForResults=False
            if not waitingForResults:
                # adapt to the next best level: the one with the most relative change
                currentElectrons = self.adaptiveGeneratorElectrons.getCurrentResult()
                currentIons = self.adaptiveGeneratorIons.getCurrentResult()
                
                # identify the next most relevant levels
                levelElectrons = self.adaptiveGeneratorElectrons.getMostRelevant()
                levelIons = self.adaptiveGeneratorIons.getMostRelevant()
                deltaElectrons = self.adaptiveGeneratorElectrons.getDelta(levelElectrons)
                deltaIons = self.adaptiveGeneratorIons.getDelta(levelIons)                
                
                #print(list(i for i in levelElectrons),list(i for i in levelIons))
                #print(self.adaptiveGeneratorElectrons.getRelevanceOfActiveSet())
                #print(self.adaptiveGeneratorIons.getRelevanceOfActiveSet())

                # compare relative differences of the different levels
                if abs(deltaElectrons / currentElectrons) > abs(deltaIons / currentIons):
                    # choose which one based on the relative change to the current value
                    l_adapt=levelElectrons
                    adaptedByElectrons=True
                    print("adapted by electrons " + str(deltaElectrons) + " vs " + str(deltaIons))
                else:
                    l_adapt=levelIons
                    adaptedByElectrons=False
                    print("adapted by ions " + str(deltaIons) + " vs " + str(deltaElectrons))

                # adapt this level for all generators
                self.adaptiveGeneratorElectrons.adaptLevel(l_adapt)
                self.adaptiveGeneratorIons.adaptLevel(l_adapt)
                self.adaptiveSEMElectrons.adaptLevel(l_adapt)
                self.adaptiveSEMIons.adaptLevel(l_adapt)

                # print output so we see what happened
                try:
                    print("adapted (" + thingToString(l_adapt) + "): [" \
                          + str(self.adaptiveGeneratorElectrons.getCurrentResult()) \
                          + "; "  + str(self.adaptiveGeneratorIons.getCurrentResult()) \
                          + "] +- [" + str(self.adaptiveSEMElectrons.getCurrentResult()) \
                          + "; "  + str(self.adaptiveSEMIons.getCurrentResult()) + "]"\
                    )
                    #print(len(self.adaptiveGeneratorElectrons.getOldSet()), \
                    #      len(self.adaptiveGeneratorElectrons.getRelevanceOfActiveSet()))
        
                    # update totalSEM and absAdaptedDelta to check the termination criterion
                    if adaptedByElectrons:
                        totalSEM=self.adaptiveSEMElectrons.getCurrentResult()
                        absAdaptedDelta=abs(self.adaptiveGeneratorElectrons.getDelta(levelElectrons))
                    else:
                        totalSEM=self.adaptiveSEMIons.getCurrentResult()
                        absAdaptedDelta=abs(self.adaptiveGeneratorIons.getDelta(levelIons))
                    print("total SEM " + str(totalSEM) + " adapted delta " + str(absAdaptedDelta)) 
                except RuntimeError:
                    # means that smaller-level results are missing
                    pass
                    
        print(len(self.adaptiveGeneratorElectrons.getOldSet()) < numGrids, (totalSEM < absAdaptedDelta), \
                not waitingForResults, waitingForNumberOfResults < 5, \
                not(len(self.adaptiveGeneratorElectrons.getRelevanceOfActiveSet()) == 0))
        if(len(self.adaptiveGeneratorElectrons.getOldSet()) > numGrids):
            print("Terminated because desired number of grids was reached ("+str(numGrids)+"): we are done.")
        if (totalSEM > self.termination_criterion_sem_factor * absAdaptedDelta):
            print("Terminated because stopping criterion was reached (SEM of "+str(totalSEM)+\
                    " is higher than delta of "+str(absAdaptedDelta)+\
                    "): maybe the adaptation continues when simulations have run longer -- check your job queue.")
        if (waitingForResults or waitingForNumberOfResults > 4 or (len(self.adaptiveGeneratorElectrons.getRelevanceOfActiveSet()) == 0)):
            print("Terminated because waiting for results of simulation runs: we are not done yet.")

        # print output 
        try:
            print("final result: " + str(self.adaptiveGeneratorElectrons.getCurrentResult()) \
                          + "; "  + str(self.adaptiveGeneratorIons.getCurrentResult())\
                          + "] +- [" + str(self.adaptiveSEMElectrons.getCurrentResult()) \
                          + "; "  + str(self.adaptiveSEMIons.getCurrentResult()) + "]"\
            )
        except RuntimeError:
            # means that smaller-level results are missing
            pass
        
        oldSet = self.adaptiveGeneratorElectrons.getOldSet()
        print("adapted levels: " + thingToString(oldSet))
        self.adaptiveGeneratorElectrons.adaptAllKnown()
        activeAndOldSet = self.adaptiveGeneratorElectrons.getOldSet()

        oldCoefficients = pysgpp.getStandardCoefficientsFromLevelSet(oldSet)
        allCoefficients = pysgpp.getStandardCoefficientsFromLevelSet(activeAndOldSet)
        
        schemeToCSV(thingToList(oldSet), thingToList(oldCoefficients), 'oldSet.csv')
        schemeToCSV(thingToList(activeAndOldSet), thingToList(allCoefficients), 'allSet.csv')
        
        if (self.output3d):
            try:
                # cf. https://stackoverflow.com/questions/1482308/how-to-get-all-subsets-of-a-set-powerset
                from itertools import chain, combinations
                s = ['l_0', 'l_1', 'l_2', 'l_3', 'l_4']
                combinations = list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))
                combinations = [c for c in list(combinations) if len(c) == 3]
                for c in combinations:
                    pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveGeneratorElectrons, c, "e_"+str(c)+".pdf")
                    pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveGeneratorIons, c, "i_"+str(c)+".pdf")
                    pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveSEMElectrons, c, "ee_"+str(c)+".pdf")
                    pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveSEMIons, c, "ei_"+str(c)+".pdf")
            except Exception as e:
                print(e)
                pass


def schemeToCSV(levels, coefficients, name):
    dflist = []
    for level, coefficient in zip(levels, coefficients):
        d = {}
        for dim in range(len(levels[0])):
            key = str('l_'+str(dim))
            d[key] = level[dim]
        d['probname'] = sim_launcher.l_vec_to_string(level)
        d['coefficient'] = coefficient
        dataframe = pd.DataFrame(d, index=[0])
        dflist.append(dataframe)
    dataframe = pd.concat(dflist, ignore_index=True)
    dataframe.to_csv(name)

