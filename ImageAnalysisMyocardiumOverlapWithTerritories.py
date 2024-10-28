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
		if self.Args.OutputFolder is None:
		#	if self.InputFileName.find("/")>=0:self.Args.OutputFolder=self.Args.InputFileName(self.Args.InputFileName.replace("/")[-1],"")
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


		#Read the Labels and Extract Territories with given label on command line
		print ("------ Reading the Territory Labels")
		if self.Args.InputFileName.find("/")>=0:infile=open(self.Args.InputFileName.replace(self.Args.InputFileName.split("/")[-1],"MBF_Territories_Labels.dat"),'r')
		else: infile=open("MBF_Territories_Labels.dat",'r')
		infile.readline()
		TerritoryLabels=[]
		for LINE in infile:
			line=LINE.split()
			if line[1].find(self.Args.TerritoryTag)>=0: TerritoryLabels.append(int(line[0]))
		infile.close()

		#Create a new set of array
		ThresholdArray=np.zeros(Volume.GetNumberOfPoints())
		for i in range(Volume.GetNumberOfPoints()):
			if int(Volume.GetPointData().GetArray("TerritoryMaps").GetValue(i)) in TerritoryLabels:
				ThresholdArray[i]=1
		
		#Store the data in Surface array
		#Tags for out or inner surface
		ThresholdArrayVTK=numpy_to_vtk(ThresholdArray,deep=True)
		ThresholdArrayVTK.SetName("TerritoryLabels_%s"%self.Args.TerritoryTag)
		Volume.GetPointData().AddArray(ThresholdArrayVTK)
		Volume.Modified()

		#Apply the Threshold to Obtain Territory Region
		TerritoryVolume=ThresholdByUpper(Volume,"TerritoryLabels_%s"%self.Args.TerritoryTag,1)
		WriteVTUFile(self.Args.OutputFolder+"/MBF_Territory_%s.vtu"%self.Args.TerritoryTag,TerritoryVolume)

		#Apply the Constant Threshold to Obtain Volume
		ThresholdVolume=ThresholdInBetween(Volume,self.Args.ArrayName,0,self.Args.LowerThreshold)
		
		#Loop over the TerritoryVolume to See how many points are not in Threshold Volume
		ThresholdVolumeCoords=vtk_to_numpy(ThresholdVolume.GetPoints().GetData())
		counter=0
		for i in range(TerritoryVolume.GetNumberOfPoints()):
			if TerritoryVolume.GetPoint(i) in ThresholdVolumeCoords:
				counter+=1
		print (counter)
			
		
		




if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will objectively extract the ischemic lesion by using FWHM from core of ischemic lesion to outer regions with normal perfusion.")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the myocardial volume.")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName", default="ImageScalars", help="The array name where the data is stored")
	
	parser.add_argument('-LowerThreshold', '--LowerThreshold', type=float, required=True, dest="LowerThreshold", default="LowerThreshold", help="The lower threshold to identify ischemia.")
	
	parser.add_argument('-TerritoryTag', '--TerritoryTag', type=str, required=True, dest="TerritoryTag", default="TerritoryTag", help="The territory that will be .")
	
	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")

	args=parser.parse_args()
	ImageAnalysisMyocardiumOverlapWithTerritories(args).Main()
