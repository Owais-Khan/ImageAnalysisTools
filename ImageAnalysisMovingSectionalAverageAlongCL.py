#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:09:54 2023

@author: ana
"""


import vtk
import argparse
from utilities import ReadVTPFile, vtk_to_numpy, ReadVTUFile, numpy_to_vtk, WriteVTPFile, WriteVTUFile, ThresholdInBetween
import numpy as np
import os
import re


class ImageAnalysisMovingAverageAlongCL():
    def __init__(self, Args):
        self.Args = Args
        self.volume = ReadVTUFile(self.Args.InputVolumeFile)
        self.num = re.findall(r'\d+', self.Args.InputVolumeFile)
        

    #defining a clipping function
    def SphereClipp(self, center, Radius, input, p):
        #define the Sphere
        Sphere = vtk.vtkSphere()
        Sphere.SetCenter(center)
        Sphere.SetRadius(Radius*0.8)
      

        #Implement vtkclipping filter "sphere"
        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(input)
        clipper.SetClipFunction(Sphere)
        clipper.InsideOutOn()
        #clipper.GetOutputInformation(1)
        clipper.Update()
        #WriteVTUFile(f"./clipper{self.num[-1]}/clipper{self.num[-1]}_{p}.vtu", clipper.GetOutput())
        #WriteVTUFile("clipper.vtu", clipper.GetOutput())
      
        #AverageFilter on each Sphere
        SphereOutput = clipper.GetOutput().GetPointData().GetArray("scalars")
        SphereOutput = vtk_to_numpy(SphereOutput)
        averagePixelValue = np.average(SphereOutput)
        return averagePixelValue

    
    def Main(self):
        
        #Read the centerlinefile and calculate the Average MaximumInscribedSphereRadius
        CLfile = ReadVTPFile(self.Args.InputCLFile)
        SphereRadiusArray = CLfile.GetPointData().GetArray("MaximumInscribedSphereRadius")
        SphereRadiusArray = vtk_to_numpy(SphereRadiusArray)
        Npoints = CLfile.GetNumberOfPoints()
        PixelValueAvg = np.zeros((Npoints-2,1))
        os.system(f"vmtkcenterlinemeshsections -centerlinesfile {self.Args.InputCLFile} -ifile {self.Args.InputVolumeFile} -ofile abc.vtp")
        SurfaceSection = ReadVTPFile("./abc.vtp")
        #os.system(f"mkdir ./clipper{self.num[-1]}")
        for p in range(Npoints-2):
            center = CLfile.GetPoint(Npoints-(p+1))
            Radius = SphereRadiusArray[Npoints-(p+1)]

            #Clip a Sphere out of the Volume
            #PixelValueAvg[p] = ImageAnalysisMovingAverageAlongCL(args).SphereClipp(center,Radius,self.volume)
            
            #Clip a Circle out of the Section
            section_ = ThresholdInBetween(SurfaceSection, "SectionIds", p, p)
            print(f"{p+1}/{Npoints-2}")
            #SectionOutput = section_.GetPointData().GetArray("scalars")
            #SectionOutput = vtk_to_numpy(SectionOutput)
            PixelValueAvg[p] = ImageAnalysisMovingAverageAlongCL(args).SphereClipp(center,Radius,section_,p)
            
        #Project the AveragePixelValue to the centerline and a text file
        PixelValueAvg = np.flipud(PixelValueAvg)
        PixelValueAvg = np.append(PixelValueAvg,0)
        PixelValueAvg = np.append(PixelValueAvg,0)
        PixelValueAvgVTK = numpy_to_vtk(PixelValueAvg)
        PixelValueAvgVTK.SetName("AvgPixelValue")
        CLfile.GetPointData().AddArray(PixelValueAvgVTK)
        OutputVolumeFile = f'OutputFolder/CenterLine_{self.num[-1]}.vtp'
        WriteVTPFile(OutputVolumeFile, CLfile)
        os.system("rm -rf ./abc.vtp")
            
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the volume file having pixel value projected on and the center-line of the lumen and retearns the average pixel values")
    
    #InputFiles
    parser.add_argument('-InputVolumeFile', '--InputVolumeFile', type = str, required = True, dest = "InputVolumeFile", help = "the volume file having pixel value projected on")
    parser.add_argument('-InputCLFile', '--InputCLFile', type = str, required = True, dest = "InputCLFile", help = "the center-line of the lumen")
    #parser.add_argument('-InputSurfaceFile', '--InputSurfaceFile', type = str, required = True, dest = "InputSurfaceFile", help = "the surface file")
    
    args = parser.parse_args()
    ImageAnalysisMovingAverageAlongCL(args).Main()