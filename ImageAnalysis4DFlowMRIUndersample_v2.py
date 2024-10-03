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
			InputFolder=self.Args.InputFolder.split("/")[-1]
			self.Args.OutputFolder=self.Args.InputFolder.replace(InputFolder,InputFolder+"_Undersampled")
			os.system("mkdir %s"%self.Args.OutputFolder)
			print ("---The OutputFolder is: %s"%self.Args.OutputFolder)

	def Main(self):
		#Get the number of time steps
		InputFiles=sorted(glob((self.Args.InputFolder+"/*.vtu")))
		N_ts=len(InputFiles)
		if N_ts is None: 
			print ("No .vtu files found in the input folder. Exiting...")
			exit(1)

		#Load the first file in the folder and get image parameters
		print (InputFiles)
		Data_=ReadVTIFile(InputFiles[0])
		N_pts=Data_.GetNumberOfPoints()
		Spacing   =Data_.GetSpacing()
		SpacingMag=(Spacing[0]**2+Spacing[1]**2+Spacing[2]**2)**0.5
		Origin    =Data_.GetOrigin()
		Dims=Data_.GetDimensions()
               
		#Create VTK Image Data to store undersampled data
		DataUS=vtk.vtkImageData()
		DataUS.SetDimensions((int(Dims[0]/self.Args.Factor),int(Dims[1]/self.Args.Factor), int(Dims[2]/self.Args.Factor)))
		DataUS.SetSpacing((Spacing[0]*self.Args.Factor,Spacing[1]*self.Args.Factor,Spacing[2]*self.Args.Factor))
		DataUS.SetOrigin(Origin)
		DataUS.AllocateScalars(vtk.VTK_FLOAT,1) #Use 2 for complex image
		DataUS.Modified() 
       
		#Get an Array of U, V and W to Store Data
		U_US=np.zeros(shape=(N_pts,N_ts))
		V_US=np.zeros(shape=(N_pts,N_ts))
		W_US=np.zeros(shape=(N_pts,N_ts))

		#Create a cube of under-sample dimensions
		Sphere=CreateSphere((0,0,0),SpacingMag*self.Args.Factor)

                

		#Loop over undersample data and average values both spatially and temporally
		for TS in range(0,N_ts):
			Data_=ReadVTIFile(InputFiles[TS])
			for i in range(0,DataUS.GetNumberOfPoints()):
				#Create a Sphere
				Sphere.SetCenter(DataUS.GetPoint(i))
				#Apply Clipper to the High-Resolution Dataset	
				Clipper.SetInputData(Data_)
				Clipper.Update()
				sphere_=Clipper.GetOutput()
				WriteVTUFile("abc.vtu",sphere_)
				print (i)	
				if i==5000:

					exit(1)
			

	
			
			
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will spatio-temporally undersample the data by a given factor (E.g., see Fatahi et al., CMPB, 2020)")
        
	parser.add_argument('-InputFolder', '--InputFolder', required=True, dest="InputFolder",help="A folder that contains VTK or VTI files.")
	
	parser.add_argument('-Factor', '--Factor', required=False, dest="Factor", default=2, help="Undersampling factor in space and time.")
        
	parser.add_argument('-OutputFolder', '--OutputFolder', required=False, dest="OutputFolder",help="The folder name to store the undersample data.")
                        
	args=parser.parse_args()
	ImageAnalysis4DFlowMRIUndersample(args).Main()

