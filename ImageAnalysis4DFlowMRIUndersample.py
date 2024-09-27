import os
import numpy as np
import scipy.io
import vtk
import argparse
from glob import glob
from utilities import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class ImageAnalysis4DFlowMRIUndersample():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			InputFolder=self.Args.InputFileName.split("/")[-1]
			self.Args.OutputFolder=self.Args.InputFolder.replace(InputFolder,"_Undersampled")
			os.system("mkdir %s"%self.Args.OutputFolder)
			print ("---The OutputFolder is: %s"%self.Args.OutputFolder)

	def Main(self):
		#Get the number of time steps
		InputFiles=glob(sorted(self.Args.InputFolder+"/*.vti"))
		N_ts=len(InputFiles)
		if N_ts is None: 
			print ("No .vti files found in the input folder. Exiting...")
			exit(1)

		#Load the first file in the folder
		Data_=ReadVTIFile(InputFiles[0])
		print (dir(Data_)) 
	
			
			
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will spatio-temporally undersample the data by a given factor.")
        
	parser.add_argument('-InputFileName', '--InputFolder', required=True, dest="InputFolder",help="A folder that contains VTK or VTI files.")
	
	parser.add_argument('-Factor', '--Factor', required=False, dest="Factor", default=2, help="Undersampling factor in space and time.")
        
	parser.add_argument('-OutputFolder', '--OutputFolder', required=False, dest="OutputFolder",help="The folder name to store the undersample data.")
                        
	args=parser.parse_args()
	ImageAnalysis4DFlowMRIUndersample(args).Main()

