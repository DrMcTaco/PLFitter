# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:06:21 2016

@author: Dan
"""
spectrumFolder = 'test_data' #Name of the folder with the spectra in it.
#Needs to be in the same directory as this script.

import lorrentzianFitter as lf
import os


with open(spectrumFolder + '\\FitParams.txt', 'w') as fHandle:#Setup Fit Parameter output file
    fHandle.write('FileName\tA-Intensity\tA-Width\tA-PeakEnergy\tAIntensity\tAWidth\tAPeakEnergy\tBIntensity\tBWidth\tBPeakEnergy\tLinearM\tLinearB\tA-Max\tAMax\tBMax\tAbsoluteMaxIntensity\tAbsoluteMaxEnergy\n')

with open(spectrumFolder + '\\FitParams.txt', 'ab') as fHandle:
    with open(spectrumFolder + '\\refit.txt', 'a') as fHandle2:
        for file in os.listdir(spectrumFolder):#do for each file in Spectra folder
            print('Fitting ' + file, end ='')
            try: 
                lf.PLFit(spectrumFolder, file, fHandle)
                print(' done!')
            except RuntimeError as err:
                print('.\n Runtime error. Likely time out.')
                fHandle2.write(file + '\n')
            except lf.FileTypeError as err:
                print('. Not a data file. Skipped.')
                