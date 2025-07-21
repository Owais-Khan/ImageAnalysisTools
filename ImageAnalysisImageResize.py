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
		InputFiles=sorted(glob((self.Args.InputFolder+"/*.vti")))
		N_ts=len(InputFiles)
		if N_ts is None: 
			print ("No .vtk files found in the input folder. Exiting...")
			exit(1)

		#Load the first file in the folder and get image parameters
		Data_=ReadVTIFile(InputFiles[0])
		
		#Resize the Image
		resize = vtk.vtkImageResize()
		resize.SetInputData(Data_)
		resize.SetResizeMethodToMagnificationFactors()
		resize.SetMagnificationFactors(self.Args.Factor,self.Args.Factor,self.Args.Factor)
		resize.InterpolateOn()
		resize.Update()
		DataResized_=resize.GetOutput()
	
		#Create Velocity Array
		VelocityResized_=np.zeros(shape=(DataResized_.GetNumberOfPoints(),3))
		
	
		#Loop over all the points and extract velocity
		ClipRadius=np.linalg.norm(DataResized_.GetSpacing())/2.
		
		#Loop over all the points
		for j in range(DataResized_.GetNumberOfPoints()):
			print (j,DataResized_.GetNumberOfPoints())
			Point__=DataResized_.GetPoint(j)
			Sphere__=CreateSphere(Point__,ClipRadius)
	
			#Clip a sphere	
			ClipData__=vtk.vtkClipDataSet()
			ClipData__.SetInputData(Data_)
			ClipData__.SetClipFunction(Sphere__)
			ClipData__.InsideOutOn()
			ClipData__.Update()
			ClipData__=ClipData__.GetOutput()	
			#Average Velocity within the Clipped Section			
			Velocity__=vtk_to_numpy(ClipData__.GetPointData().GetArray("Velocity"))
			VelocityX__=np.average(Velocity__[:,0])
			VelocityY__=np.average(Velocity__[:,1])
			VelocityZ__=np.average(Velocity__[:,2])
			
			VelocityResized_[j]=[VelocityX__,VelocityY__,VelocityZ__]

		VelocityResizedVTK_=numpy_to_vtk(VelocityResized_)	

		WriteVTUFile(os.path.join(self.Args.OutputFolder,"ClipData_.vtu"),ClipData__.GetOutput())
	
			
			
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will spatio-temporally undersample the data by a given factor (E.g., see Fatahi et al., CMPB, 2020)")
        
	parser.add_argument('-InputFolder', '--InputFolder', required=True, dest="InputFolder",help="A folder that contains VTK or VTI files.")
	
	parser.add_argument('-Factor', '--Factor', required=True, dest="Factor", type=float, default=2, help="Undersampling factor in space and time.")
        
	parser.add_argument('-OutputFolder', '--OutputFolder', required=False, dest="OutputFolder",help="The folder name to store the undersample data.")
                        
	args=parser.parse_args()
	ImageAnalysis4DFlowMRIUndersample(args).Main()

