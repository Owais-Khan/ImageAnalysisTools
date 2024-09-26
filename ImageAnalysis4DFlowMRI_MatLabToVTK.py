import os
import numpy as np
import scipy.io
import vtk
import argparse
from glob import glob
from utilities import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class ImageAnalysis4DFlowMRI_MatLabToVTK():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			InputFileName=self.Args.InputFileName.split("/")[-1]
			self.Args.OutputFolder=self.Args.InputFileName.replace(InputFileName,"")+"./VTK_Files"
			os.system("mkdir %s"%self.Args.OutputFolder)
			print ("---The OutputFolder is: %s"%self.Args.OutputFolder)

	def Main(self):
		print ("Reading the Matlab File...")
		Data = scipy.io.loadmat(self.Args.InputFileName)
		PCMR = np.array(Data["PCMR"])
		
		#Number of Timesteps
		No_TimeSteps=int(max(PCMR[:,3]))
		print ("Number of Timesteps: %d"%No_TimeSteps)
		
		#Number of Coordinates
		No_Coordinates=int(float(len(Data["PCMR"]))/No_TimeSteps)
		print ("Number of Coordinates: %d"%No_Coordinates)
		
		#Create a separate split Array
		PCMR=np.array_split(PCMR,No_TimeSteps)

		#Spatial Resolution
		SpacingX=float(Data["spatialres"][0].split("x")[0])	
		SpacingY=float(Data["spatialres"][0].split("x")[0])	
		SpacingZ=float(Data["spatialres"][0].split("x")[0])
		print ("Space X Y Z is: %.03f %.03f %.03f"%(SpacingX, SpacingY, SpacingZ))		


		#Dimensions
		DimX=int((max(PCMR[0][:,0])-min(PCMR[0][:,0]))/SpacingX)+1
		DimY=int((max(PCMR[0][:,1])-min(PCMR[0][:,1]))/SpacingY)+1
		DimZ=int((max(PCMR[0][:,2])-min(PCMR[0][:,2]))/SpacingZ)+1
		print ("Dimensions in X Y and Z is: %d %d %d"%(DimX,DimY,DimZ))	

		#Get Matching IDs between Matlab Array and VTK Image
		CoordsMatLab=np.zeros(shape=(No_Coordinates,3))
		for i in range(No_Coordinates): 
			CoordsMatLab[i,0]=PCMR[0][i][0]
			CoordsMatLab[i,1]=PCMR[0][i][1]
			CoordsMatLab[i,2]=PCMR[0][i][2]

		print ("Min and Max X from MatLab is:",min(CoordsMatLab[:,0]),max(CoordsMatLab[:,0]))
		print ("Min and Max Y from MatLab is:",min(CoordsMatLab[:,1]),max(CoordsMatLab[:,1]))
		print ("Min and Max Z from MatLab is:",min(CoordsMatLab[:,2]),max(CoordsMatLab[:,2]))

		#Create VTK Image Data
		ImageData=vtk.vtkImageData()
		ImageData.SetDimensions(DimX,DimZ,DimY)
		ImageData.SetSpacing(SpacingX,SpacingY,SpacingZ)
		ImageData.SetOrigin(min(CoordsMatLab[:,0]),min(CoordsMatLab[:,2]),min(CoordsMatLab[:,1]))
		ImageData.AllocateScalars(vtk.VTK_FLOAT,1) #Use 2 for complex image
		
		#Get Nearest Point IDs
		No_Coords_VTK=ImageData.GetNumberOfPoints()
		if No_Coords_VTK!=No_Coordinates:
			print ("The number of MatLab Coordinates do not match VTK Image Coordinates")
			print ("Exiting...")
		CoordsVTK=np.zeros(shape=(No_Coordinates,3))
		for i in range(No_Coords_VTK):
			CoordsVTK[i,:]=ImageData.GetPoint(i)[:]	

		print ("Min and Max X from VTK is:",min(CoordsVTK[:,0]),max(CoordsVTK[:,0]))
		print ("Min and Max Y from VTK is:",min(CoordsVTK[:,1]),max(CoordsVTK[:,1]))
		print ("Min and Max Z from VTK is:",min(CoordsVTK[:,2]),max(CoordsVTK[:,2]))


	
		for i in range(0,No_TimeSteps):
			print ("------ Looping over: %d"%i)
			velocity_=np.zeros(shape=(No_Coordinates,3))
			counter=0
			for j in range(No_Coordinates):		
				velocity_[j,0]=PCMR[i][counter,-3]	
				velocity_[j,1]=PCMR[i][counter,-2]	
				velocity_[j,2]=PCMR[i][counter,-1]
				counter+=1	
				velmag_=(velocity_[j,0]**2+velocity_[j,1]**2+velocity_[j,2]**2)**0.5
				ImageData.GetPointData().GetArray("ImageScalars").SetValue(j,velmag_)	
				 
			velocityVTK_=numpy_to_vtk(velocity_)
			velocityVTK_.SetName("Velocity")
			ImageData.GetPointData().AddArray(velocityVTK_)
			ImageData.Modified()
			WriteVTIFile(self.Args.OutputFolder+"/Velocity_%.03d.vti"%i,ImageData)
		
			
			
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take in a volumetric dataset with X Y Z T Mag U V W array written in matlab format.")
        
	parser.add_argument('-InputFileName', '--InputFileName', required=True, dest="InputFileName",help="A Matlab file.")
        
	parser.add_argument('-OutputFolder', '--OutputFolder', required=False, dest="OutputFolder",help="The folder name where all the output velocities will be stored.")
                        
	args=parser.parse_args()
	ImageAnalysis4DFlowMRI_MatLabToVTK(args).Main()

