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

class ImageAnalysisMyocardiumStatisticsIschemicRegion():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=self.Args.InputFolder
		if self.Args.InputFileName is None:
			self.Args.InputFileName=os.path.join(self.Args.InputFolder,"MBF_Territories.vtu")
			print ("Tesselation File: %s"%self.Args.InputFileName)
		if self.Args.InputFileName2 is None:
			self.Args.InputFileName2=os.path.join(self.Args.InputFolder,"MBF_IschemicTerritory_1.vtu")
			print ("Ischemic File: %s"%self.Args.InputFileName2)

	def Main(self):
		#Read the Volume
		Volume=self.ReadFile()

		#Find the array name with work "scalars"
		if self.Args.ArrayName is None:
			for i in range(Volume.GetPointData().GetNumberOfArrays()):
				ArrayName_=str(Volume.GetPointData().GetArrayName(i))
				if ArrayName_.find("calars")>0: self.Args.ArrayName=ArrayName_
		print ("--- The array that contains MBF values is: %s"%self.Args.ArrayName)

		#Read the Labels and Extract Territories with given label on command line
		print ("------ Reading the Territory Labels")
		infile=open(os.path.join(self.Args.InputFolder,"MBF_Territories_Labels.dat"),'r')
		infile.readline()
		TerritoryLabels=[]
		for LINE in infile:
			line=LINE.split()
			if line[1].find(self.Args.TerritoryTag)>=0: TerritoryLabels.append(int(line[0]))
		infile.close()

		#Create a new set of array
		print ("--- Extracting the %s Territory"%self.Args.TerritoryTag)
		ThresholdArray=np.zeros(Volume.GetNumberOfPoints())
		for i in range(Volume.GetNumberOfPoints()):
			if int(Volume.GetPointData().GetArray("TerritoryMaps").GetValue(i)) in TerritoryLabels:
				ThresholdArray[i]=1
		
		#Tag the Terrotiry that is supposed to be ischemic
		ThresholdArrayVTK=numpy_to_vtk(ThresholdArray,deep=True)
		ThresholdArrayVTK.SetName("TerritoryLabels_%s"%self.Args.TerritoryTag)
		Volume.GetPointData().AddArray(ThresholdArrayVTK)
		Volume.Modified()

		#Apply the Threshold to Obtain Territory Region
		print ("--- Saving the %s Territory"%self.Args.TerritoryTag)
		TerritoryTesselation=ThresholdByUpper(Volume,"TerritoryLabels_%s"%self.Args.TerritoryTag,1)
		WriteVTUFile(self.Args.OutputFolder+"/MBF_Territory_%s.vtu"%self.Args.TerritoryTag,TerritoryTesselation)
		#----------------------------------------------------------------

		#Compute the Statistics for True Ischemic and Tesselation-detected Ischemic Territory		
		#Read the True Ischemic Territory
		TerritoryIschemic=ReadVTUFile(self.Args.InputFileName2)
		TerritoryIschemic=LargestConnectedRegion(TerritoryIschemic) #Extract Largest Connected Region
		TerritoryIschemicSurface=ExtractSurface(TerritoryIschemic)
		
		
		print ("--- Computing the Statistics...")
		StatisticsEntireLV=Statistics(Volume,self.Args.ArrayName)	
		StatisticsIschemic=Statistics(TerritoryIschemic,self.Args.ArrayName,NormalizationValue=StatisticsEntireLV["75thPerct"])
		StatisticsTesselation=Statistics(TerritoryTesselation,self.Args.ArrayName,NormalizationValue=StatisticsEntireLV["75thPerct"])



		#Check for overlap points
		TerritoryTesselationPoints=vtk_to_numpy(TerritoryTesselation.GetPoints().GetData())
		TerritoryIschemicPoints=vtk_to_numpy(TerritoryIschemic.GetPoints().GetData())

		print ("--- Computing Overlap")	
		#Convert into string
		TerritoryTesselationString=[]
		TerritoryIschemicString=[]
		for i in range(TerritoryTesselation.GetNumberOfPoints()): TerritoryTesselationString.append(str(TerritoryTesselationPoints[i]))
		for i in range(TerritoryIschemic.GetNumberOfPoints()): TerritoryIschemicString.append(str(TerritoryIschemicPoints[i]))
		TotalPointsString=TerritoryTesselationString+TerritoryIschemicString
		UniquePointsString=set(TotalPointsString)
		OverlapCounter=len(TotalPointsString)-len(UniquePointsString)


                #Calculate the DICE score for the overlap
		print ("--- Points in Ischemic Territory:    %d"%TerritoryIschemic.GetNumberOfPoints())
		print ("--- Points in Tesselation Territory: %d"%TerritoryTesselation.GetNumberOfPoints())	
		print ("--- Overlap Points:                  %d"%OverlapCounter)	
	
		#Compute the Dice Score
		StatisticsTesselation["DiceScore"]=(2*OverlapCounter)/(TerritoryTesselation.GetNumberOfPoints()+TerritoryIschemic.GetNumberOfPoints())
	
		#Write Output to the File
		OutputFileName=self.Args.InputFileName2.replace(".vtu","_IschemicOverlapStatistics.dat")
		print ("--- Writing the Statistics: %s"%OutputFileName)
		outfile=open(OutputFileName,'w')
		outfile.write("Entire LV Statistics\n")
		for key in StatisticsEntireLV:   
			if StatisticsEntireLV[key] is not None:outfile.write ("%s: %.05f\n"%(key,StatisticsEntireLV[key]))
			else: outfile.write("%s: None\n"%key)
		outfile.write("-"*25+"\n")
		outfile.write("Tesselation Statistics\n")
		for key in StatisticsTesselation:  
			if StatisticsTesselation[key] is not None: outfile.write ("%s: %.05f\n"%(key,StatisticsTesselation[key]))
			else: outfile.write("%s: None\n"%key)
		outfile.write("-"*25+"\n")
		
		outfile.write("Threshold Statistics\n")
		for key in StatisticsIschemic: 
			if StatisticsIschemic[key] is not None: outfile.write ("%s: %.05f\n"%(key,StatisticsIschemic[key]))
			else: outfile.write("%s: None\n"%key)	
		outfile.write("-"*25+"\n")
		outfile.close()
        
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
	parser = argparse.ArgumentParser(description="This script compare the true ischemic territory with tessellation-detected ischemic territory")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The input folder where all of the files are stroed.")
	
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=False, dest="InputFileName",help="The vtu file that contains the myocardial volume with territory labels.")
	
	parser.add_argument('-InputFileName2', '--InputFileName2', type=str, required=False, dest="InputFileName2",help="The vtu file that contains the true ischemic territory.")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName", help="The array name where the data is stored")
	
	parser.add_argument('-TerritoryTag', '--TerritoryTag', type=str, required=True, dest="TerritoryTag", default="TerritoryTag", help="The territory that will be .")
	
	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")

	args=parser.parse_args()
	ImageAnalysisMyocardiumStatisticsIschemicRegion(args).Main()
