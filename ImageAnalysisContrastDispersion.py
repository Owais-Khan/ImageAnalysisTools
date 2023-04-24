#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 20:33:17 2023

@author: ana
"""


import numpy as np
import glob
import argparse
from math import sqrt, atan, pi
from utilities import ReadVTPFile
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

class ImageAnalysisAdvectionDiffusionAlongCL():
    def __init__(self,Args):
        self.Args = Args
        
    def MAFilter(self):
        input_ = self.x_denoised
        window_length = 20
        FilteredSignal = np.zeros(np.size(input_)-window_length)
        for point in np.size(input_):
            FilteredSignal[point] = np.average(input_[point:point+window_length])
        return FilteredSignal
        
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
        PixelValArray = PixelValArray[:-3,0:7] 
        self.PixelValArray = PixelValArray
        
        #Calculating Velocity
        #D = 0.01
        point_1 = CLFile.GetPoint(1)
        point_2 = CLFile.GetPoint(2)
        x_step = sqrt((point_1[0] - point_2[0])**2+(point_1[1] - point_2[1])**2+(point_1[2] - point_2[2])**2) #Calculate the distance between the CL points
        
        print(x_step)
        time_step = 2.69#4 #assuming the heart-rate to be 60-per-min and the pictures are taken every 4-cycle
        (nx, nt) = PixelValArray.shape
        x_arr = np.arange(0, nx*x_step, x_step)
        t_arr = np.arange(0, nt*time_step, time_step)
        dc_dx = 0
        dc_dt = 0
        for x_ in range(nx):
            model = LinearRegression()
            model.fit(t_arr.reshape(-1, 1),PixelValArray[x_,:].reshape(-1, 1))
            dc_dt += model.coef_[0][0]
            '''plt.scatter(t_arr, PixelValArray[x_,:], color = 'black')
            pred = model.predict(t_arr.reshape(-1, 1))
            plt.plot(t_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
            plt.text(6, 250, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 41)
            plt.title("Contrast over time")
            plt.xlabel("time (sec)")
            plt.ylabel("Contrast (HU)")
            plt.show()
            exit(1)'''
        dc_dt /= nx
        x_denoised = PixelValArray.sum(axis = 1)/nx
        self.x_denoised = x_denoised
        x_MAFiltered = ImageAnalysisAdvectionDiffusionAlongCL().MAFilter()
        model = LinearRegression()
        model.fit(x_arr.reshape(-1, 1),x_MAFiltered.reshape(-1, 1))
        dc_dx += model.coef_[0][0]
        plt.scatter(x_arr, x_MAFiltered, color = 'black')
        pred = model.predict(x_arr.reshape(-1, 1))
        plt.plot(x_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
        plt.text(12, 183, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 2)
        plt.title("Averaged Contrast over Space")
        plt.xlabel("distance from inflow cap (mm)")
        plt.ylabel("Contrast (HU)")
        plt.axis([0, 30, 150, 200])
        plt.show()
        '''
        for t_ in range(nt):
            model = LinearRegression()
            model.fit(x_arr.reshape(-1, 1),PixelValArray[:,t_].reshape(-1, 1))
            dc_dx += model.coef_[0][0]
            plt.scatter(x_arr, PixelValArray[:,t_], color = 'black')
            pred = model.predict(x_arr.reshape(-1, 1))
            plt.plot(x_arr.reshape(-1, 1), pred.reshape(-1, 1), color = 'red')
            plt.text(12, 54.8, f"slope = {int(model.coef_[0][0]*100)/100}", rotation = 24)
            plt.title("Contrast over space")
            plt.xlabel("distance from inflow cap (mm)")
            plt.ylabel("Contrast (HU)")
            plt.show()'''
        dc_dx /= nt
        velocity = dc_dt/dc_dx
        print(abs(velocity), 'mm/sec')
        
        
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the centerline files with the averaged pixel value array and calculates the advection diffusion velocity")
    
    #InputFiles
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = "InputFolder", help = "The Folder containing the centerline files with averaged pixel values")
    
    args = parser.parse_args()
    ImageAnalysisAdvectionDiffusionAlongCL(args).Main()
    