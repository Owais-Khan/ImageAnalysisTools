import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
import argparse
from utilities import *
from time import time
class ImageAnalysisLabelImage():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Load the VTI Image
		print ("--- Loading the Image: %s"%self.Args.InputFileName)
		Image=ReadVTIFile(self.Args.InputFileName)
	
		#Load the Surface file contain a closed model of the aorta-coronary vasculature
		print ("--- Loading the Surface: %s"%self.Args.InputSurface)
		Surface=ReadVTPFile(self.Args.InputSurface)


		#Creata a points array
		print ("--- Converting Image into Point Array. This may take some time...")
		PointsVTK=vtk.vtkPoints()
		PointsVTK.SetNumberOfPoints(Image.GetNumberOfPoints())
		progress_=0
		for i in range(Image.GetNumberOfPoints()):
			PointsVTK.SetPoint(i,Image.GetPoint(i))
			progress_=PrintProgress(i,Image.GetNumberOfPoints(),progress_)
		
		print ("--- Converting Image Points into a Polydata")
		#Convert into a polydata format
		pdata_points = vtk.vtkPolyData()
		pdata_points.SetPoints(PointsVTK)
	
		print ("--- Checking whether the Image Points are inside the Surface")	
		#Loop over all of the points in the image and check if they fall inside the closed surface.
		selectEnclosed = vtk.vtkSelectEnclosedPoints()
		selectEnclosed.SetInputData(pdata_points) #Points in the Image
		selectEnclosed.SetSurfaceData(Surface) #Surface Model
		selectEnclosed.SetTolerance(0.00000000001)
		selectEnclosed.Update()

		print ("--- Copying Image Labels to a Numpy Array")
		#Yes or No for all of the points in the Image
		ImageLabels=vtk_to_numpy(selectEnclosed.GetOutput().GetPointData().GetArray('SelectedPoints'))

		#ImageLabels=np.zeros(Image.GetNumberOfPoints())
		#for i in range(Image.GetNumberOfPoints()):
		#	ImageLabels[i]=selectEnclosed.GetOutput().GetPointData().GetArray('SelectedPoints').GetTuple(i)[0]

		print ("--- Appending Numpy Image Labels to the VTK Image Volume")
		print ("------ 1=aorta, 2=left coronary tree and 3== right coronary tree and 0=none") 
		ImageNew=SurfaceAddArray(Image,ImageLabels,"Labels")

		print ("--- Writing Image to the Output file: %s"%self.Args.OutputFileName)
		WriteVTIFile(self.Args.OutputFileName,ImageNew)	


if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will label the coronary CTA based on 1=aorta, 2=left coronary tree and 3=coronary tree, and 0=none coronary/aorta vasculature.")

        #Input filename of the coronary CTA
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="File name of the vti image containing the coronary CTA")

	#Input filename of the coronary segmented surface.
	parser.add_argument('-InputSurface', '--InputSurface', type=str, required=True, dest="InputSurface",help="File name of the input surface generated from simvascular, including the caps. This is the mesh-complete.exterior.vtp file")
	
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=True, dest="OutputFileName",help="File name in which to store the Image along with the Labels array.")
    
	args=parser.parse_args()
	ImageAnalysisLabelImage(args).Main()

