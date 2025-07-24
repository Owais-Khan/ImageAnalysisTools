import os
import numpy as np
import vtk
import argparse
from glob import glob
from utilities import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

import subprocess
class ImageAnalysisPointCloudToSurface():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=os.path.join(self.Args.InputFolder,"NiftiLabels")
			subprocess.run("mkdir %s"%self.Args.OutputFolder, shell=True)
	
		self.Args.Tags=[str(value) for value in self.Args.Tags.split(",")] #separate the tags
		print ("\n"+"-"*30)
		print ("The FileName to Assign Masks are:",self.Args.Tags)

	def Main(self):
		#Read the Image File
		print ("\n"+"-"*30)
		print ("---Reading the SimVascular Input Folder: %s"%self.Args.InputFolder)
		ImageFiles=sorted(glob(os.path.join(self.Args.InputFolder,"Images/*.vti")))
		ImageFile=ImageFiles[0]
		if len(ImageFiles)==0: raise Exception("No .vti files found in the image folder.")
		
		#Output filename
		self.Args.OutputFileName=os.path.join(self.Args.OutputFolder,os.path.basename(ImageFile))
		print ("------ Reading Image: %s"%ImageFiles[0])
		ImageVTI=ReadVTIFile(ImageFiles[0])

		# Create an intersection filter
		mesh_complete_filename=sorted(glob(os.path.join(self.Args.InputFolder,"Meshes/*mesh-complete/mesh-complete.mesh.vtu")))[0]
		print ("------ Reading mesh-complete.mesh.vtu file: %s"%mesh_complete_filename)
		MeshCompleteVTP=ReadVTUFile(mesh_complete_filename)

			
		print("\n"+"-"*30)
		#Get the Tags and Coordinates based on mesh-complete folder
		Coords,Tags=self.FindTags() #Contains tags (1,2,3 etc) and Coordinates. 

		#Build Static Cell Locator of the Surface 
		CellLocator=vtk.vtkStaticCellLocator()
		CellLocator.SetDataSet(MeshCompleteVTP)
		CellLocator.BuildLocator()	
	
		#Loop over the volume and get points
		NumberOfImagePoints=ImageVTI.GetNumberOfPoints()

		#Make the Scalar_ Array 0 to make background 0
		ScalarsVTK=numpy_to_vtk(np.zeros(NumberOfImagePoints),deep=True)
		ScalarsVTK.SetName("Scalars_")
		ImageVTI.GetPointData().AddArray(ScalarsVTK)
	
		print ("\n --- Looping over the Image volume to assign tags. This may take some time...")	
		#Loop over the entire volume and fill the tags according to vessel name
		for i in range (NumberOfImagePoints):
			point_=ImageVTI.GetPoint(i)
			test_=CellLocator.FindCell(point_)
			if test_ != -1: #If image point is found inside the mesh
				value_,id_,min_dist_=ClosestPoint(point_,Coords)
				ImageVTI.GetPointData().GetArray("Scalars_").SetValue(i,Tags[id_])
							
		print ("------ Writing File: %s.nii.gz"%self.Args.OutputFileName)
		WriteNiftiFile(self.Args.OutputFileName.replace(".vti",".nii.gz"),ImageVTI)



	def FindTags(self):
		#Read all of the surface files
		print ("\n--- Reading Surface Files from Mesh Folder")
		SurfaceFiles=sorted(glob(os.path.join(self.Args.InputFolder,"Meshes/*mesh-complete/mesh-surfaces/wall_*.vtp")))
		if len(SurfaceFiles)==0: raise Exception("No Surfaces found in the mesh-complete folder")
		else: print ("------Number of Surfaces found: %d"%len(SurfaceFiles))

		#Loop over all of the surfaces
		Coords=[]
		Tags=[]
		currentTag=0 #0=Aorta, 1=Left Coronary, 2=Right Coronary
		for SurfaceFile_ in SurfaceFiles:
			#Read the surface files
			print ("------ Looping over: %s"%os.path.basename(SurfaceFile_))
			FileData_=ReadVTPFile(SurfaceFile_)
			FileName_=os.path.basename(SurfaceFile_)

			#Loop over the tag names
			for i in range(len(self.Args.Tags)): 
				if FileName_.find(self.Args.Tags[i])>=0: currentTag=i+1 #Assign 0=background, 1=tag1, 2=tag2 
			print ("--------- Label Id is: %d"%currentTag)
			
			coords_=vtk_to_numpy(FileData_.GetPoints().GetData())
			for value in coords_: 
				Coords.append(value)
				Tags.append(currentTag)
		Coords=np.array(Coords)
		Tags=np.array(Tags)

			
		return Coords,Tags
            
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will convert SimVascular folder to Nifti format. ")
	parser.add_argument('-InputFolder', '--InputFolder', required=True, dest="InputFolder",help="The mesh complete folder.")
	parser.add_argument("-Tags", "--Tags", required=False, default="+", help="The list of arguments that contain labels.")
	parser.add_argument('-OutputFolder', '--OutputFolder', required=False, dest="OutputFolder",help="Folder to store the labelled maps.")
                        
	args=parser.parse_args()
	ImageAnalysisPointCloudToSurface(args).Main()

