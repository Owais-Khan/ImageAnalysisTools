#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:32:29 2023

@author: ana
"""

import numpy as np
from math import pi, sqrt
from scipy.fft import fft2, fftfreq, ifft2
import glob
import argparse
from utilities import ReadVTPFile
import pandas as pd

class ImageAnalysisAdvectionDiffusionAlongCL():
    def __init__(self,Args):
        self.Args = Args
        
    def Main(self):
        
        #Store the cl file names inside the input folder
        FolderName = self.Args.InputFolder
        FileNames  = glob.glob(f'{FolderName}/*.vtp')
        FileNames  = sorted(FileNames)
        NFiles = FileNames.__len__()
        Npoints = ReadVTPFile(FileNames[0]).GetNumberOfPoints()
        
        #Read the averaged pixel values of each cl file and storing them into an excel sheet
        PixelValArray = np.zeros((Npoints, NFiles))
        for i in range(NFiles):
            CLFile = ReadVTPFile(FileNames[i])
            PixelValArray[:,i] = CLFile.GetPointData().GetArray("AvgPixelValue")
            
        PixelValArray = np.delete(PixelValArray, (-1), axis = 0)
        PixelValArray = np.delete(PixelValArray, (-1), axis = 0)
        df = pd.DataFrame(PixelValArray)
        df.to_excel(excel_writer = f"{self.Args.InputFolder}/AveragePixelValues.xlsx")
        
        #Calculating Velocity
        D = 0.01
        point_1 = CLFile.GetPoint(1)
        point_2 = CLFile.GetPoint(2)
        x_step = sqrt((point_1[0] - point_2[0])**2+(point_1[1] - point_2[1])**2+(point_1[2] - point_2[2])**2) #Calculate the distance between the CL points
        #print(x_step)
        SphereVolume = 4/3*pi*(0.8*x_step)**3
        time_step = 4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        Contrast_fft = np.array(fft2(PixelValArray))
        (nx, nt) = PixelValArray.shape
        (fx, ft) = (fftfreq(nx, d=SphereVolume), fftfreq(nt, d=time_step))
        dC_dt_fft = 1j*2*pi*ft*Contrast_fft
        dC_dx_fft = 1j*2*pi*fx*Contrast_fft.transpose()
        d2C_dx2_fft = 1j*2*pi*fx*dC_dx_fft
        dC_dx_fft = dC_dx_fft.transpose()
        d2C_dx2_fft = d2C_dx2_fft.transpose()
        
        #velocity_fft = np.dot((D*d2C_dx2_fft - dC_dt_fft),dC_dx_fft)
        velocity_fft = (D*d2C_dx2_fft - dC_dt_fft)/dC_dx_fft

        #Back to the time-space domain: Velocity Along the Centerline
        velocity = ifft2(velocity_fft[1:,:])
        
        df1 = pd.DataFrame(abs(velocity))
        df1.to_excel(excel_writer = f"{self.Args.InputFolder}/VelocityAlongCL.xlsx")
        
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the centerline files with the averaged pixel value array and calculates the advection diffusion velocity")
    
    #InputFiles
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = "InputFolder", help = "The Folder containing the centerline files with averaged pixel values")
    
    args = parser.parse_args()
    ImageAnalysisAdvectionDiffusionAlongCL(args).Main()
    
    