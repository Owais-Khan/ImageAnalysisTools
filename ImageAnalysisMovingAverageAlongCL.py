#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 11:47:54 2023

@author: ana
"""

import vtk
import argparse
from utilities import ReadVTPFile, vtk_to_numpy, ReadVTUFile, numpy_to_vtk, WriteVTPFile
import numpy as np
import os
import re

class ImageAnalysisMovingAverageAlongCL():
    def __init__(self, Args):
        self.Args = Args
        self.volume = ReadVTUFile(self.Args.InputVolumeFile)
        self.num = re.findall(r'\d+', self.Args.InputVolumeFile)
        

    #defining a clipping function
    def SphereClipp(self, center, Radius):
        #define the Sphere
        Sphere = vtk.vtkSphere()
        Sphere.SetCenter(center)
        Sphere.SetRadius(Radius*0.8)

            
        #Implement vtkclipping filter "sphere"
        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(self.volume)
        clipper.SetClipFunction(Sphere)
        clipper.InsideOutOn()
        clipper.GetOutputInformation(1)
        clipper.Update()
            
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
        SphereRadiusAverage = np.average(SphereRadiusArray)
        
        #Extract the centerline of the lumen with a new resamplingstep
        CL_File_Name = self.Args.InputVolumeFile.replace(".vtu","_cl.vtp")
        os.system(f"vmtkcenterlines -ifile {self.Args.InputSurfaceFile} -ofile {CL_File_Name} -endpoints 1 -resampling 1 -resamplingstep {SphereRadiusAverage}")
        
        #looping over the center-line points
        CLfile = ReadVTPFile(CL_File_Name)
        SphereRadiusArray = CLfile.GetPointData().GetArray("MaximumInscribedSphereRadius")
        SphereRadiusArray = vtk_to_numpy(SphereRadiusArray)
        Npoints = CLfile.GetNumberOfPoints()
        PixelValueAvg = np.zeros((Npoints,1))
        for p in range(Npoints):
            center = CLfile.GetPoint(p)
            Radius = SphereRadiusArray[p]
            PixelValueAvg[p] = ImageAnalysisMovingAverageAlongCL(args).SphereClipp(center,Radius)
           
        #Project the AveragePixelValue to the centerline and a text file
        PixelValueAvgVTK = numpy_to_vtk(PixelValueAvg)
        PixelValueAvgVTK.SetName("AvgPixelValue")
        CLfile.GetPointData().AddArray(PixelValueAvgVTK)
        OutputVolumeFile = f'OutputFolder/CenterLine{self.num[-1]}.vtp'
        WriteVTPFile(OutputVolumeFile, CLfile)
        
            
if __name__ == "__main__":
    #Description
    parser = argparse.ArgumentParser(description = "This script gets the volume file having pixel value projected on and the center-line of the lumen and retearns the average pixel values")
    
    #InputFiles
    parser.add_argument('-InputVolumeFile', '--InputVolumeFile', type = str, required = True, dest = "InputVolumeFile", help = "the volume file having pixel value projected on")
    parser.add_argument('-InputCLFile', '--InputCLFile', type = str, required = True, dest = "InputCLFile", help = "the center-line of the lumen")
    parser.add_argument('-InputSurfaceFile', '--InputSurfaceFile', type = str, required = True, dest = "InputSurfaceFile", help = "the surface file")
    
    args = parser.parse_args()
    ImageAnalysisMovingAverageAlongCL(args).Main()