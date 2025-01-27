import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *
import csv
from scipy import stats

class ImageAnalysisSimVascularCoronaryMaps():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
                #Read the vtu file
		print ("--- Reading Mesh: %s"%self.Args.SimVasc)
		Volume=self.ReadFile()

                #Find the array name with word "scalars"
		if self.Args.ArrayName is None:
			for i in range(VentricleMesh.GetPointData().GetNumberOfArrays()):
				ArrayName_=str(VentricleMesh.GetPointData().GetArrayName(i))
				if ArrayName_.find("calars")>0: self.Args.ArrayName=ArrayName_


		#Apply the Constant Threshold to Obtain Volume
		ThresholdVolumeIschemic=ThresholdInBetween(Volume,self.Args.ArrayName,0,self.Args.LowerThreshold)
		ThresholdVolumeHealthy=ThresholdInBetween(Volume,self.Args.ArrayName,self.Args.LowerThreshold,10000000)

		#Now extract the surface to find the closest distance
		print ("--- Extracting the surface")
		SurfaceIschemic=ExtractSurface(ThresholdVolumeIschemic)
		SurfaceHealthy=ExtractSurface(ThresholdVolumeHealthy)

		#Compute the Distance Between Ischemic Surface and Entire Myocardium
		DistanceIschemic=vtk.vtkImplicitPolyDataDistance()
		DistanceIschemic.SetInput(SurfaceIschemic)
		DistanceHealthy =vtk.vtkImplicitPolyDataDistance()
		DistanceHealthy.SetInput(SurfaceHealthy)


		#Create an array to store ischemic vs non-ischemic tags
		IschemicMyocardiumTag=np.zeros(Volume.GetNumberOfPoints())
		NumberOfPoints=Volume.GetNumberOfPoints()
		progress_old_=-1
		print ("--- Assigning the Ischemic and Healthy Territories")
		for i in range(Volume.GetNumberOfPoints()):
			#Print the progress
			progress_=PrintProgress(i,NumberOfPoints,progress_old_)
			progress_old_=progress_
			#Find the Distance closest to the Surface	
			dist_ischemic_=DistanceIschemic.EvaluateFunction(Volume.GetPoint(i))
			dist_healthy_ =DistanceHealthy.EvaluateFunction(Volume.GetPoint(i))
			if dist_ischemic_<=dist_healthy_:
				IschemicMyocardiumTag[i]=1
		 
		IschemicMyocardiumTagVTK=numpy_to_vtk(IschemicMyocardiumTag,deep=True)
		IschemicMyocardiumTagVTK.SetName("IschemicTerritory")
		Volume.GetPointData().AddArray(IschemicMyocardiumTagVTK)
		Volume.Modified()


		print ("--- Threshold Ischemic Volume with Epicardial/Endocardial Walls")
		#Apply Threshold to Obtain Ischemic Region
		ThresholdVolumeIschemic=ThresholdInBetween(Volume,"IschemicTerritory",1,1)
		
		#Find the largest connected volume
		print ("--- Extracting the Largest Connected Region")
		ConnectedVolume=vtk.vtkConnectivityFilter()
		ConnectedVolume.SetInputData(ThresholdVolumeIschemic)
		ConnectedVolume.SetExtractionModeToLargestRegion()
		ConnectedVolume.Update()
		ConnectedVolumeData=ConnectedVolume.GetOutput()

		OutputFileName=self.Args.InputFileName.replace(".vtu","_%d_IschemicTerritory.vtu"%self.Args.LowerThreshold)
		print ("--- Writing the Output File: %s"%OutputFileName)
		WriteVTUFile(OutputFileName,ConnectedVolumeData)	
		
		
	def ReadFile(self):
		print ("--- Reading %s"%self.Args.InputFileName)
		if   self.Args.InputFileName[-4:]==".vtk":
			Volume=ReadVTKFile(self.Args.InputFileName)   #Read a VTK volume stack
			Volume=ThresholdByUpper(Volume,self.Args.ArrayName,1) #Convert to an unstructured grid 
		elif self.Args.InputFileName[-4:]==".vtu":
			Volume=ReadVTUFile(self.Args.InputFileName) #Read a VTU unstructured grid
		else:
			print ("The extension %s is not valid for volume"%self.Args.InputFileName[-4:])
			print ("Exiting...")
			exit(1)
		return Volume


if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take a SimVascular Folder and create image.nii.gz and label.nii.gz files for nnUNet trainning. This script is specifically designed for coronary arteries. The labelled maps will be 0=background, 1=aorta, 2=LCA, 3=RCA.")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFileName",help="The volumetric mesh in the mesh-complete folder.")
	
	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")

	args=parser.parse_args()
	ImageAnalysisMyocardiumOverlapWithTerritories(args).Main()
