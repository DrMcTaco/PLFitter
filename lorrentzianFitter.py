# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 16:19:51 2016

@author: Dan Rubin
"""

import numpy as np
import spc
import pyqtgraph as pg
from scipy.optimize import curve_fit
import pyqtgraph.exporters

def lorrentz1(x, a, b, c):#Single Lorrentzian function
    return (a*b/(2.*np.pi))/((x-c)**2+(b/2.)**2)
       
def lorrentz3(x, a, b, c, r, s, t, h, k, l): #Sum of threee lorrentzian functions
    return lorrentz1(x, a, b, c) + lorrentz1(x, r, s, t) + lorrentz1(x, h, k, l)

def linear(x, m, b):
    return m*x + b

def fitFunc(x, a, b, c, r, s, t, h, k, l, m, n):
    return lorrentz3(x, a, b, c, r, s, t, h, k, l) + linear(x, m, n)

def loadData(folder, fileName):
    if fileName.endswith(".spc") or fileName.endswith('.SPC'):
        f = spc.File(folder + '\\' + fileName)
        ev = f.x
        counts = f.sub[0].y
    elif fileName.endswith(".txt") or fileName.endswith(".csv") or fileName.endswith(".PRN") or fileName.endswith(".prn"):
        ev, counts = np.genfromtxt(folder + '\\' + fileName, skip_footer=1, unpack=True, usecols =(0,1))  #load data from text file
    else:
        raise Exception
    
    return ev, counts

def PLFit(folder, fileName, bestFitFile):    
    eV, counts = loadData(folder, fileName)
    
    

    if eV[0] > eV[1]:#Ensure the eV array is in ascending order. reverse order if incorrect
        eV = np.flipud(eV)
        counts = np.flipud(counts)
    
    limits = ([0., 0., 1.79, 0., 0., 1.79, 0., 0., 1.89, -np.inf, -np.inf],[np.inf, 1, 1.9, np.inf, 1, 1.9, np.inf, .5, 2.05, np.inf, np.inf])#bounds on best fit parameters
    bestFitParameters, covariance = curve_fit(fitFunc, eV, counts, p0=(1, .1, 1.84, 1, .1, 1.88, 1, .01, 2.02, 500, -500), bounds = limits, method='trf', max_nfev=10000) #fit data to triple lorrentzian

    A = []#initialize output arrays
    B = []
    A_prime = []
    Spectrum_fit = []
    line = []

    for e in eV: # build output arrays
        A.append(lorrentz1(e, bestFitParameters[0], bestFitParameters[1], bestFitParameters[2]))
        B.append(lorrentz1(e, bestFitParameters[3], bestFitParameters[4], bestFitParameters[5]))
        A_prime.append(lorrentz1(e, bestFitParameters[6], bestFitParameters[7], bestFitParameters[8]))
        line.append(linear(e, bestFitParameters[9], bestFitParameters[10]))
        Spectrum_fit.append(fitFunc(e, bestFitParameters[0], bestFitParameters[1], bestFitParameters[2], bestFitParameters[3], bestFitParameters[4], bestFitParameters[5], bestFitParameters[6], bestFitParameters[7], bestFitParameters[8],bestFitParameters[9],bestFitParameters[10]))

    bestFitParameters.resize(len(bestFitParameters)+5)
    bestFitParameters[11] = (bestFitParameters[0]*bestFitParameters[1])/(2*np.pi*(bestFitParameters[1]/2)**2)
    bestFitParameters[12] = (bestFitParameters[3]*bestFitParameters[4])/(2*np.pi*(bestFitParameters[4]/2)**2)
    bestFitParameters[13] = (bestFitParameters[6]*bestFitParameters[7])/(2*np.pi*(bestFitParameters[7]/2)**2)
    bestFitParameters[14] = max(Spectrum_fit)
    bestFitParameters[15] = eV[Spectrum_fit.index(max(Spectrum_fit))]
    
# Why all this nonsense is here: https://stackoverflow.com/questions/21203907/numpy-savetxt-heterogenous-data
    name = np.array(fileName)
    dt = np.dtype([('name', name.dtype.str)] + [('data', bestFitParameters.dtype.str, bestFitParameters.shape[0])])
    temp = np.empty(1, dtype=dt)
    temp['data'] = bestFitParameters.tolist()
    temp['name'] =  fileName
    dtall = np.dtype([('name', name.dtype.str)] + bestFitParameters.dtype.descr*bestFitParameters.shape[0])
    temp2 = temp.view(dtall)      
    id_fmt = '%s'   # The format of the id.
    float_fmt = '%.6f'  # The floating point format.
    fmt = id_fmt + bestFitParameters.shape[0] * (' ' + float_fmt)    
    np.savetxt(bestFitFile, temp2, delimiter = '\t', fmt=fmt)
    
    bestFitParameters.resize(len(eV))
    
    if fileName.endswith(".spc") or fileName.endswith('.SPC'):
        np.savetxt(folder + '\\' + fileName.split('.')[0] + '.txt', np.transpose([eV, counts, Spectrum_fit, A, B, A_prime, line, bestFitParameters]), delimiter = '\t', fmt='%.6f')# save to text file
    else:
        np.savetxt(folder + '\\' + fileName, np.transpose([eV, counts, Spectrum_fit, A, B, A_prime, line, bestFitParameters]), delimiter = '\t', fmt='%.6f')# save to text file
    
    data = [counts, Spectrum_fit, A, B, A_prime, line]
    plotData(eV, data, fileName, folder)
    
    
def plotData(x, data, fileName, folder):
    pw = pg.plot()
    pw.setTitle(fileName)

    for i in range(len(data)):
        pw.plot(x, data[i], pen = (i ,len(data)))
        
    exp = pg.exporters.ImageExporter(pw.plotItem)
    exp.parameters()['width'] = 640
    exp.export(folder + '\\' + fileName.split('.')[0] + '.png')


 

if __name__ == "__main__":   
    with open('Files_Peak_Analysis' + '\\FitParams.txt', 'w') as fHandle:
        fHandle.write('FileName\tA-Intensity\tA-Width\tA-PeakEnergy\tAIntensity\tAWidth\tAPeakEnergy\tBIntensity\tBWidth\tBPeakEnergy\tLinearM\tLinearB\tA-Max\tAMax\tBMax\tAbsoluteMaxIntensity\tAbsoluteMaxEnergy\n')    
    with open('Files_Peak_Analysis' + '\\FitParams.txt', 'ab') as fHandle:
        PLFit('Files_Peak_Analysis', 'test.spc', fHandle) #for testing purpose.
