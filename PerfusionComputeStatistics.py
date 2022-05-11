import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *
import matplotlib.pyplot as plt
import csv

class PerfusionComputeStatistics():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			os.system("mkdir Results")
			self.Args.OutputFolder="./Results/"

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

		
		#Compute All of the Scalar quantities
		print ("--- Reading Scalar Values")
		Npts=Volume.GetNumberOfPoints()
		DataArray=np.zeros(Npts)
		for i in range(Npts):
			DataArray[i]=Volume.GetPointData().GetArray(self.Args.ArrayName).GetValue(i)
			
		#Perform Statistics
		Mean=np.average(DataArray)
		Stdev=np.average(DataArray)
		Median50=np.percentile(DataArray,50)
		Median75=np.percentile(DataArray,75)
		Median80=np.percentile(DataArray,80)
		Median85=np.percentile(DataArray,85)
		Median90=np.percentile(DataArray,90)
		Median95=np.percentile(DataArray,95)

		print ("Writing Statistics: %s/MBF_Statistics.dat"%self.Args.OutputFolder)
		#Write the data to a file
		outfile=open("%s/MBF_Statistics.dat"%self.Args.OutputFolder,'w')
		outfile.write("Mean    MBF: %.05f\n"%Mean)
		outfile.write("Stdev   MBF: %.05f\n"%Stdev)
		outfile.write("Median  MBF: %.05f\n"%Median50)
		outfile.write("Perc75  MBF: %.05f\n"%Median75)
		outfile.write("Perc80  MBF: %.05f\n"%Median80)
		outfile.write("Perc85  MBF: %.05f\n"%Median85)
		outfile.write("Perc90  MBF: %.05f\n"%Median90)
		outfile.write("Perc95  MBF: %.05f\n"%Median95)

		outfile.close()

		print ("Writing Raw Data: %s/MBF_Raw.dat"%self.Args.OutputFolder)
		outfile=open("%s/MBF_Raw.txt"%self.Args.OutputFolder,'w')
		for i in range(Npts):
			outfile.write("%d %.05f\n"%(i,DataArray[i]))
		outfile.close()

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute data statistics on MBF volume")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the myocardial volume")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName", default="ImageScalars", help="The array name where the data is stored")

	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")

	args=parser.parse_args()
	PerfusionComputeStatistics(args).Main()
