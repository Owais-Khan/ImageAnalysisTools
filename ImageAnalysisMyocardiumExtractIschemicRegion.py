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

class ImageAnalysisMyocardiumOverlapWithTerritories():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.InputFileName.find("/")>=0:
			FileName=self.Args.InputFileName.split("/")[-1]
			self.Args.OutputFolder=self.Args.InputFileName.replace(FileName,"")
		else:	
			self.Args.OutputFolder="./"

	def Main(self):
                #Read the vtu file
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

		#Apply the Constant Threshold to Obtain Volume
		ThresholdVolumeIschemic=ThresholdInBetween(Volume,self.Args.ArrayName,0,self.Args.LowerThreshold)
		ThresholdVolumeNormal=ThresholdInBetween(Volume,self.Args.ArrayName,self.Args.LowerThreshold,10000000)

		#Extract the largest connect region
		print ("--- Extracting the largest connected region of the ischemic myocardium.")
		ConnectivityFilterIschemic=vtk.vtkConnectivityFilter()
		ConnectivityFilterIschemic.SetInputData(ThresholdVolumeIschemic)
		ConnectivityFilterIschemic.SetExtractionModeToLargestRegion()
		ConnectivityFilterIschemic.Update()
		ConnectivityFilterIschemic=ConnectivityFilterIschemic.GetOutput()

		print ("--- Extracting the largest connected region of the normal myocardium.")
		ConnectivityFilterNormal=vtk.vtkConnectivityFilter()
		ConnectivityFilterNormal.SetInputData(ThresholdVolumeNormal)
		ConnectivityFilterNormal.SetExtractionModeToLargestRegion()
		ConnectivityFilterNormal.Update()
		ConnectivityFilterNormal=ConnectivityFilterNormal.GetOutput()

		#Now extract the surface to find the closest distance
		print ("--- Extracting the surface")
		SurfaceIschemic=ExtractSurface(ConnectivityFilterIschemic)
		SurfaceNormal=ExtractSurface(ConnectivityFilterNormal)

		#Separate all of the coordinates of the volume either into ischemic or normal myocardium
		print ("--- Separating the Myocardium into Ischemic and Normal Territories")
		IschemicMyocardiumTag=np.zeros(Volume.GetNumberOfPoints())
		SurfaceIschemicCoords=vtk_to_numpy(SurfaceIschemic.GetPoints().GetData())
		SurfaceNormalCoords  =vtk_to_numpy(SurfaceNormal.GetPoints().GetData())

		for i in range(Volume.GetNumberOfPoints()):
			coord_=Volume.GetPoint(i)
			value_,argument_,min_dist_Ischemic=ClosestPoint(coord_,SurfaceIschemicCoords)
			value_,argument_,min_dist_Normal=ClosestPoint(coord_,SurfaceNormalCoords)
			if min_dist_Ischemic<=min_dist_Normal:
				IschemicMyocardiumTag[i]=1


		 
		IschemicMyocardiumTagVTK=numpy_to_vtk(IschemicMyocardiumTag,deep=True)
		IschemicMyocardiumTagVTK.SetName("IschemicTerritory")
		Volume.GetPointData().AddArray(IschemicMyocardiumTagVTK)
		Volume.Modified()

		print ("--- Writing the Output File: %s"%self.Args.OutputFolder)
		WriteVTUFile(self.Args.OutputFolder+"/abc2.vtu",Volume)	
		
		




if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will objectively extract a single volume of ischemic region using threshold.")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the myocardial volume.")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName", default="ImageScalars", help="The array name where the data is stored")
	
	parser.add_argument('-LowerThreshold', '--LowerThreshold', type=float, required=True, dest="LowerThreshold", default="LowerThreshold", help="The lower threshold to identify ischemia.")
	
	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")

	args=parser.parse_args()
	ImageAnalysisMyocardiumOverlapWithTerritories(args).Main()
