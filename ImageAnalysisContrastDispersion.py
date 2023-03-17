#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 20:33:17 2023

@author: ana
"""


import numpy as np
import glob
import argparse
from math import sqrt
from utilities import ReadVTPFile
from sklearn.linear_model import LinearRegression

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
            
        #RadiusArray = CLFile.GetPointData().GetArray("MaximumInscribedSphereRadius")
        '''
        PixelValArray = np.delete(PixelValArray, (-1), axis = 0)
        PixelValArray = np.delete(PixelValArray, (-1), axis = 0)
        PixelValArray = np.delete(PixelValArray, (-1), axis = 0)
        '''
        # Taking the linear part of the lumen and the upslope samples
        PixelValArray = PixelValArray[11:34,6:9] 
        
        
        #Calculating Velocity
        #D = 0.01
        point_1 = CLFile.GetPoint(1)
        point_2 = CLFile.GetPoint(2)
        x_step = sqrt((point_1[0] - point_2[0])**2+(point_1[1] - point_2[1])**2+(point_1[2] - point_2[2])**2) #Calculate the distance between the CL points
        #print(x_step)
        time_step = 4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (x, t) = PixelValArray.shape
        x_arr = np.arange(0, x*x_step, x_step)
        t_arr = np.arange(0, t*time_step, time_step)
        dc_dx = 0
        dc_dt = 0
        for x_ in range(x):
            model = LinearRegression()
            model.fit(t_arr.reshape(-1, 1),PixelValArray[x_,:].reshape(-1, 1))
            dc_dt += model.coef_
        dc_dt /= x
        for t_ in range(t):
            model = LinearRegression()
            model.fit(x_arr.reshape(-1, 1),PixelValArray[:,t_].reshape(-1, 1))
            dc_dx += model.coef_
        dc_dx /= t
        velocity = dc_dt/dc_dx
        print(velocity)
        
        
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the centerline files with the averaged pixel value array and calculates the advection diffusion velocity")
    
    #InputFiles
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = "InputFolder", help = "The Folder containing the centerline files with averaged pixel values")
    
    args = parser.parse_args()
    ImageAnalysisAdvectionDiffusionAlongCL(args).Main()
    