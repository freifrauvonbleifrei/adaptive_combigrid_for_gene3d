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


def get_results_string(result, sem):
    if len(result) == 1:
        return "[" + str(result[0])\
                    + "] +- [" + str(sem[0]) + "]"\

    if len(result) == 2:
        return "[" + str(result[0]) \
                    + "; "  + str(result[1])\
                    + "] +- [" + str(sem[0]) \
                    + "; "  + str(sem[1]) + "]"\

class DimensionalAdaptation:
    def __init__(self, lmin, lmax, omega=1., \
                prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'),\
                gene_path=os.environ.get('ADAPTATION_GENE3D_PATH'),\
                output3d=False, skiplevels=[]):
        self.prob_prepath=prob_prepath
        self.gene_path=gene_path
        self.output3d=output3d
        self.lmin=lmin
        self.lmax=lmax
        self.skiplevels=skiplevels
        self.num_species=Qes_data.get_num_species()
        self.minimum_time_length=float(os.environ.get('ADAPTATION_PARAM_MIN_RUNTIME'))
        #self.minimum_time_length=500
        self.termination_criterion_sem_factor=float(os.environ.get('ADAPTATION_PARAM_TERMINATION_CRITERION_SEM_FACTOR'))
        #self.termination_criterion_sem_factor=1.5
        self.minimum_nwin=int(os.environ.get('ADAPTATION_PARAM_SEM_MINIMUM_NUMBER_WINDOWS'))
        #self.minimum_nwin=3

        dim = len(lmin)
        for skiplevel in skiplevels:
            assert(len(skiplevel) == dim)
        lmin = pysgpp.LevelVector(lmin)
        levels = pysgpp.LevelVectorVector(1,lmin)
        weightedRelevanceCalculator=pysgpp.WeightedRelevanceCalculator(omega)
        #TODO allow passing relevance generator with own omega -- for now setting omega to 1. in sgpp manually
        assert(omega == 1.)
        self.adaptiveGeneratorIons = pysgpp.AdaptiveCombinationGridGenerator(levels)#, weightedRelevanceCalculator)
        if self.num_species > 1:
            self.adaptiveGeneratorElectrons = pysgpp.AdaptiveCombinationGridGenerator(levels)#, weightedRelevanceCalculator)

        # set start results
        qes_data = Qes_data.Qes_data(sim_launcher.l_vec_longer(lmin), skiplevels=self.skiplevels)
        result = qes_data.get_result()
        sgppActiveLevelVector=pysgpp.LevelVector(lmin)
        self.adaptiveGeneratorIons.setQoIInformation(sgppActiveLevelVector, result[0])
        if self.num_species > 1:
            self.adaptiveGeneratorElectrons.setQoIInformation(sgppActiveLevelVector, result[1])

        #same for SEMs
        sem = qes_data.get_time_error()
        self.adaptiveSEMIons = pysgpp.AdaptiveCombinationGridGenerator.fromFunctionPointer(levels, \
                                        pysgpp.DoubleVector(1,sem[0]), pysgpp.squaredSummation_cb)
        delta = [self.adaptiveGeneratorIons.getDelta(sgppActiveLevelVector)]
        if self.num_species > 1:
            self.adaptiveSEMElectrons = \
                pysgpp.AdaptiveCombinationGridGenerator.fromFunctionPointer(levels, \
                                            pysgpp.DoubleVector(1,sem[1]), pysgpp.squaredSummation_cb)
            delta.append(self.adaptiveGeneratorElectrons.getDelta(sgppActiveLevelVector))
        qes_data.set_delta_to_csv(delta)

        if ((qes_data.how_long_run() < self.minimum_time_length) or qes_data.how_many_nwin_run() < self.minimum_nwin):
            #prolong simulation
            if sim_launcher.check_folder_exists(self.prob_prepath, lmin):
                if sim_launcher.check_finished(self.prob_prepath, lmin):
                    # restart the simulation
                    print("sim_launcher: restart "+str(lmin))
                    sim_launcher.restart_sim(
                        self.prob_prepath, self.gene_path, lmin)
                print("need to run the inital simulation longer!", lmin)
                raise AssertionError(
                    "need to run the inital simulation longer -- run the script again later!")
            else:
                raise AssertionError(
                    "something went wrong with the initial simulation -- check the prob_ directory")

        print("starting with: " + self.get_results_string())

    def get_results_string(self):
        if self.num_species == 1:
            return "[" + str(self.adaptiveGeneratorIons.getCurrentResult())\
                      + "] +- [" + str(self.adaptiveSEMIons.getCurrentResult()) + "]"\

        if self.num_species == 2:
            return "[" + str(self.adaptiveGeneratorIons.getCurrentResult()) \
                       + "; "  + str(self.adaptiveGeneratorElectrons.getCurrentResult())\
                      + "] +- [" + str(self.adaptiveSEMIons.getCurrentResult()) \
                      + "; "  + str(self.adaptiveSEMElectrons.getCurrentResult()) + "]"\

    def adapt_for_all(self, l_adapt):
        self.adaptiveGeneratorIons.adaptLevel(l_adapt)
        self.adaptiveSEMIons.adaptLevel(l_adapt)
        if self.num_species == 2:
            self.adaptiveGeneratorElectrons.adaptLevel(l_adapt)
            self.adaptiveSEMElectrons.adaptLevel(l_adapt)

    def run_dimadapt_algorithm(self, numGrids=10):
        absAdaptedDelta = 100000. # kinds.default_float_kind.MAX
        totalSEM=-1.
        waitingForResults=False
        waitingForNumberOfResults=0
        while (len(self.adaptiveGeneratorIons.getOldSet()) < numGrids) and \
                (totalSEM < self.termination_criterion_sem_factor* absAdaptedDelta) and \
                waitingForNumberOfResults < 10 and \
                ~(len(self.adaptiveGeneratorIons.getRelevanceOfActiveSet()) == 0):
            waitingForResults=True
            waitingForNumberOfResults=0

            # iterate the whole active set
            for activeLevelVector in self.adaptiveGeneratorIons.getActiveSet():
                #print(activeLevelVector, waitingForNumberOfResults)
                #print(not (self.adaptiveGeneratorIons.hasQoIInformation(pysgpp.LevelVector(activeLevelVector))))
                # store QoI information if we have not done it already

                if (not self.adaptiveGeneratorIons.hasQoIInformation(pysgpp.LevelVector(activeLevelVector))):
                    print(activeLevelVector, waitingForNumberOfResults)
                    qes_data = Qes_data.Qes_data(sim_launcher.l_vec_longer(activeLevelVector), skiplevels=self.skiplevels)
                    result = qes_data.get_result()
    #                result = [0.3,0.3]
                    if result:
                        waitingForResults=False
                        sgppActiveLevelVector=pysgpp.LevelVector(activeLevelVector)
                        sem=qes_data.get_time_error()
                        assert(len(sem) == self.num_species)
                        # store the results in the adaptive generators
                        self.adaptiveGeneratorIons.setQoIInformation(sgppActiveLevelVector, result[0])
                        self.adaptiveSEMIons.setQoIInformation(sgppActiveLevelVector, sem[0])
                        deltaIons = self.adaptiveGeneratorIons.getDelta(sgppActiveLevelVector)
                        delta = [deltaIons]
                        if self.num_species > 1:
                            self.adaptiveGeneratorElectrons.setQoIInformation(sgppActiveLevelVector, result[1])
                            self.adaptiveSEMElectrons.setQoIInformation(sgppActiveLevelVector, sem[1])
                            deltaElectrons = self.adaptiveGeneratorElectrons.getDelta(sgppActiveLevelVector)
                            delta = [deltaIons, deltaElectrons]
                        print("has values: " + get_results_string(result, sem))
                        if (result[0] < math.inf):
                            qes_data.set_delta_to_csv(delta)
                            if ((qes_data.how_long_run() < self.minimum_time_length) or qes_data.how_many_nwin_run() < self.minimum_nwin):
                            #prolong simulation
                                if sim_launcher.check_folder_exists(self.prob_prepath, activeLevelVector):
                                    if sim_launcher.check_finished(self.prob_prepath, activeLevelVector):
                                        # restart the simulation
                                        print("sim_launcher: restart "+str(activeLevelVector))
                                        sim_launcher.restart_sim(self.prob_prepath, self.gene_path, activeLevelVector)
                            elif any([sem[i] > 0.1 * result[i] for i in range(len(sem))]):
                            #prolong simulation
                                if sim_launcher.check_folder_exists(self.prob_prepath, activeLevelVector):
                                    if sim_launcher.check_finished(self.prob_prepath, activeLevelVector):
                                        # restart the simulation
                                        print("sim_launcher: restart because of SEM "+str(activeLevelVector))
                                        sim_launcher.restart_sim(self.prob_prepath, self.gene_path, activeLevelVector, qes_data.how_long_run() + 2 * self.minimum_time_length)

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
                # first, look at the current adaptation result
                currentIons = self.adaptiveGeneratorIons.getCurrentResult()
                if self.num_species > 1:
                    currentElectrons = self.adaptiveGeneratorElectrons.getCurrentResult()

                # identify the next most relevant levels
                # try the "relative" criterion: the delta is divided by the level's own result
                # parse shell boolean cf https://stackoverflow.com/questions/63116419/evaluate-boolean-environment-variable-in-python
                relativeAdaptationCriterion = os.environ.get('ADAPTATION_PARAM_RELATIVE_ADAPTATION', 'True').lower() in ['true', '1']
                activeLevelVectors = self.adaptiveGeneratorIons.getActiveSet()
                levelIons = self.adaptiveGeneratorIons.getMostRelevant()
                if relativeAdaptationCriterion:
                    currentAbsRelativeMaxDelta = 0.
                    for activeLevelVector in activeLevelVectors:
                        pysgppActiveLevelVector = pysgpp.LevelVector(activeLevelVector)
                        deltaActive = self.adaptiveGeneratorIons.getDelta(pysgppActiveLevelVector)
                        if not np.isnan(deltaActive):
                            relativeDeltaActive = deltaActive / \
                                              self.adaptiveGeneratorIons.getQoIInformation(pysgppActiveLevelVector)
                            if abs(relativeDeltaActive) > currentAbsRelativeMaxDelta:
                                currentAbsRelativeMaxDelta = abs(relativeDeltaActive)
                                levelIons = pysgppActiveLevelVector
                deltaIons = self.adaptiveGeneratorIons.getDelta(levelIons)

                if self.num_species > 1:
                    levelElectrons = self.adaptiveGeneratorElectrons.getMostRelevant()
                    if relativeAdaptationCriterion:
                        currentAbsRelativeMaxDelta = 0.
                        for activeLevelVector in activeLevelVectors:
                            activeLevelVector = pysgpp.LevelVector(activeLevelVector)
                            deltaActive = self.adaptiveGeneratorIons.getDelta(activeLevelVector)
                            if not np.isnan(deltaActive):
                                relativeDeltaActive = deltaActive / \
                                                   self.adaptiveGeneratorElectrons.getQoIInformation(activeLevelVector)
                                if abs(relativeDeltaActive) > currentAbsRelativeMaxDelta:
                                    currentAbsRelativeMaxDelta = abs(relativeDeltaActive)
                                    levelElectrons = activeLevelVector
                    deltaElectrons = self.adaptiveGeneratorElectrons.getDelta(levelElectrons)

                #print(list(i for i in levelElectrons),list(i for i in levelIons))
                #print(self.adaptiveGeneratorElectrons.getRelevanceOfActiveSet())
                #print(self.adaptiveGeneratorIons.getRelevanceOfActiveSet())

                # compare relative differences of the different levels, to decide which species determines the adaptation
                if self.num_species == 2 and abs(deltaElectrons / currentElectrons) > abs(deltaIons / currentIons):
                    # choose which one based on the relative change to the current value
                    l_adapt=levelElectrons
                    adaptedByElectrons=True
                    print("adapted by electrons " + str(deltaElectrons) + " vs " + str(deltaIons))
                else:
                    l_adapt=levelIons
                    adaptedByElectrons=False
                    #print("adapted by ions " + str(deltaIons) + " vs " + str(deltaElectrons))

                # adapt this level for all generators
                self.adapt_for_all(l_adapt)

                # print output so we see what happened
                try:
                    print("adapted (" + thingToString(l_adapt) + "): " \
                          + self.get_results_string() )
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

        print(len(self.adaptiveGeneratorIons.getOldSet()) < numGrids, (totalSEM < absAdaptedDelta), \
                not waitingForResults, waitingForNumberOfResults < 10, \
                not(len(self.adaptiveGeneratorIons.getRelevanceOfActiveSet()) == 0))
        if(len(self.adaptiveGeneratorIons.getOldSet()) > numGrids):
            print("Terminated because desired number of grids was reached ("+str(numGrids)+"): we are done.")
        if (totalSEM > self.termination_criterion_sem_factor * absAdaptedDelta):
            print("Terminated because stopping criterion was reached (SEM of "+str(totalSEM)+\
                    " is higher than delta of "+str(absAdaptedDelta)+\
                    "): maybe the adaptation continues when simulations have run longer -- check your job queue.")
        if (waitingForResults or waitingForNumberOfResults > 4 or (len(self.adaptiveGeneratorIons.getRelevanceOfActiveSet()) == 0)):
            print("Terminated because waiting for results of simulation runs: we are not done yet.")

        # print output 
        try:
            print("final result: " + self.get_results_string())
        except RuntimeError:
            # means that smaller-level results are missing
            pass

        oldSet = self.adaptiveGeneratorIons.getOldSet()
        print("adapted levels: " + thingToString(oldSet))
        self.adaptiveGeneratorIons.adaptAllKnown()
        activeAndOldSet = self.adaptiveGeneratorIons.getOldSet()

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
                    pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveGeneratorIons, c, "i_"+str(c)+".pdf")
                    pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveSEMIons, c, "ei_"+str(c)+".pdf")
                    if self.num_species > 1:
                        pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveGeneratorElectrons, c, "e_"+str(c)+".pdf")
                        pysgpp.extensions.combigrid.plotDeltas3d.plotDeltas3D(self.adaptiveSEMElectrons, c, "ee_"+str(c)+".pdf")
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

