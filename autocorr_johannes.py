#!/usr/bin/env python3
import os
import shutil
import math

import numpy as np
import scipy.interpolate as spi
#import matplotlib.pyplot as plt, matplotlib.cm as cm, matplotlib.colors as colors
#plt.style.use('seaborn-whitegrid')


def get_mean_error(x,y,clearDoubleTimes=True,takeCombiTimes=False,printOutput=False):
    """
    Es nimmt als input zwei 1D arrays x, y (Zeitwerte und QoI-werte) und
    spuckt die zeitgemittelte QoI aus, samt Abschätzung
    der statistischen Unsicherheit. Die Funktion lautet:
    
    Die Parameter sind:
    
    x: 1D numpy-array oder Liste mit Zeitwerten

    y: 1D numpy-array oder Liste mit QoI-werten
    
    clearDoubleTimes: Bei 'True' werden zunächst doppelt auftretende
    Zeitpunkte herausgefiltert.
    
    takeCombiTimes: Bei 'True' werden nur die doppelt auftretenden
    Zeitpunkte (Kombi-Zeitpunkte) behalten.
    
     """

    autocorrtime, autocorrindex, yeven, datamean = get_autocorr_time(x,y,clearDoubleTimes,takeCombiTimes,printOutput)

    #if printOutput:
    #    print(autocorrtime, autocorrindex, yeven, datamean)

    #calculate batch means and standard error of the means (SEM) ##TODO:make it work on irregular grid
    windowindexspan = 5*autocorrindex
    Nwin = int(len(yeven)//windowindexspan+1)
    assert Nwin > 0
    batchmeans = np.zeros(Nwin)
    for i in range(Nwin-1):
        batchmeans[i] = np.mean(yeven[i*windowindexspan:(i+1)*windowindexspan])
    if len(yeven) > (Nwin-1)*windowindexspan:
        batchmeans[Nwin-1] = np.mean(yeven[(Nwin-1)*windowindexspan:])
    else:
        Nwin -=1
        batchmeans = batchmeans[:Nwin-1]
    sbatch = np.sqrt( np.sum( (batchmeans-datamean)**2 )/Nwin/(Nwin-1) )
    if printOutput:
        print(len(yeven), str((Nwin-1)*windowindexspan))
        print(batchmeans)
        print("Nwin: "+str(Nwin))
        print("SEM: "+str(sbatch))
        print("-------------")
        print()
    return (datamean,sbatch)


def get_autocorr_time(x,y,clearDoubleTimes=True,takeCombiTimes=False,printOutput=False):
    tstart = x[0]
    tend = x[-1]
    tspan = tend - tstart

    datamean= get_mean(x,y,clearDoubleTimes,takeCombiTimes,printOutput)

    n = int(round(1*len(x)))
    xeven = np.linspace(tstart,tend,n)
    #resample to equidistant grid using (cubic) spline interpolation
    spline = spi.BSpline(*spi.splrep(x,y,k=3))
    yeven = spline(xeven)
    
    #calculate auto-correlation function
    ##g0 = np.sum((yeven-datamean)**2)
    r = np.empty_like(yeven)
    yevendiff = yeven-datamean
    r[0] = np.sum(yevendiff**2)
    for j in range(1,len(yeven)):
        r[j] = np.sum(yevendiff[j:]*yevendiff[:-j])
    r/=r[0]
    #plt.plot(xeven,r,'g.-')
    #plt.show()

    #calculate auto-correlation time
    i = 0
    thres = 1/np.exp(1)
    while r[i]>thres:
        i+=1
    autocorrindex = i
    autocorrtime = xeven[i]-tstart
    if printOutput:
        print()
        print("Auto-correlation time: "+str(autocorrtime))

    return autocorrtime, autocorrindex, yeven, datamean
    
def get_mean(x,y,clearDoubleTimes=True,takeCombiTimes=False,printOutput=False):

    # delete Points that have the same time coordinate as the next one
    if clearDoubleTimes:
        #print("Checking for double times")
        badIndices = []
        for i in range(len(x)-1):
            if x[i] == x[i+1]:
                badIndices.append(i) # append(i+1) in order to cancel the point immediately AFTER recombination
        if len(badIndices):
            newx = list(x)
            newy = list(y)
            shift = 0
            for i in badIndices:
                newx.pop(i-shift)
                newy.pop(i-shift)
                shift += 1
            x = np.array(newx)
            y = np.array(newy)
    ###############
        
    # take only Points at Recombination steps if wanted
    if takeCombiTimes:
        assert clearDoubleTimes==False,'can not take combi times while clearDoubleTimes==True'
        #print("Filtering Combi Times")
        goodIndices = []
        for i in range(len(x)-1):
            if x[i] == x[i+1]:
                goodIndices.append(i+1)
        if len(goodIndices):
            newx = []
            newy = []
            for i in goodIndices:
                newx.append(x[i])
                newy.append(y[i])
            x = np.array(newx)
            y = np.array(newy)
    ###############

    x = np.array(x)
    y = np.array(y)
    tstart = x[0]
    tend = x[-1]
    tspan = tend - tstart
    #plt.plot(x,y,'g.-')
    #plt.show()
    if printOutput:
        print("tstart = "+str(tstart))
        print("tspan = "+str(tspan))

    #resample to equidistant grid using (cubic) spline interpolation
#     spline = spi.BSpline(*spi.splrep(x,y,k=3))


    #calculate mean
    marith = np.mean(y) #normal mean
    mtrap = np.sum((y[1:]+y[:-1])*(x[1:]-x[:-1]))/(2*tspan) #quadrature using trapezoidal rule
#     mspline = spline.integrate(tstart,tend)/tspan #quadrature using spline interpolation on equidistant grid
    if printOutput:
        print()
        print("mean:         "+str(marith))
        print("mean_trap:    "+str(mtrap))
#         print("mean_spline:  "+str(mspline))
    datamean=mtrap
    return datamean

###################### EXAMPLE CODE ##########################
if __name__ == "__main__":
#    with open('/home/rentrop/Data/geneoutput/fg_nlitgaeProbs/xzseries/QesDataC.dat','w') as outfile:
#        for lx in range(7,11):
#            for lz in range(3,7):
#                nrgfile = "/home/rentrop/Data/geneoutput/fg_nlitgaeProbs/fg_nlitgae_x"+str(int(2**lx))+"z"+str(int(2**lz))+"/output/nrgC.dat"
#                if not os.path.exists(nrgfile):
#                    nrgfile = "/home/rentrop/Data/geneoutput/fg_nlitgaeProbs/fg_nlitgae_x"+str(int(2**lx))+"z"+str(int(2**lz))+"/output/nrg.dat"
#                print("File: "+str(nrgfile))
#                print()
#                if not os.path.exists(nrgfile):
#                    print("File not found!")
#                    continue
#                (mean,err) = get_mean_error(nrgfile)
#                outfile.write(str(lx)+"\t"+str(lz)+"\t"+str(mean)+"\t"+str(err)+"\n")

    nrgfile = "/lustre/cray/ws9/0/ws/ipvpolli-exahd/gene3d-dev/prob_5_5_5_5_3/out/nrg_3"
    print("File: "+str(nrgfile))
    if not os.path.exists(nrgfile):
        raise FileNotFoundError 
    x=np.arange(0,1000,0.1)
    y = np.sin(x) 
    (mean,err) = get_mean_error(x, y, printOutput=True)
