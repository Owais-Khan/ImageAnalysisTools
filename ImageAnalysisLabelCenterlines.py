import vtk
import numpy as np
import os
from glob import glob
import argparse
from utilities import *
class ImageAnalysisLabelImage():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Load the VTI Image
		print ("--- Loading the Image: %s"%self.Args.InputFileName)
		Image=ReadVTIFile(self.Args.InputFileName)
	
		#Load the centerlines
		print ("--- Loading the Centerlines: %s"%self.Args.InputSurface)
		Surface=ReadVTPFile(self.Args.InputSurface)

		#Save all of the coordinates
		CLPoints=np.array([Surface.GetPoint(i) for i in range(Surface.GetNumberOfPoints())])

		#Build a Point Locator Class
		print ("--- Creating a Cell Locator Class")
		CellLocator = vtk.vtkCellLocator()
		CellLocator.SetDataSet(Image) 
		CellLocator.BuildLocator()

		cellId = vtk.reference(0)
		c = [0.0, 0.0, 0.0]
		subId = vtk.reference(0)
		d = vtk.reference(0.0)

		#Create an empty numpy array
		CenterlineLabels=np.zeros(Image.GetNumberOfCells())
		
		#Loop over all of the points
		print ("--- Looping over centerline coordinates")
		for i in range(len(CLPoints)):
			CellId_=CellLocator.FindClosestPoint(CLPoints[i], c, cellId, subId, d)
			CenterlineLabels[cellId]=1.0
		
		#Create a new image with labels
		print ("--- Appending Numpy Image Labels of Centerline to the VTK Image Volume")
		ImageNew=SurfaceAddCellArray(Image,CenterlineLabels,"CenterlineLabels")	

		print ("--- Writing Image to the Output file: %s"%self.Args.OutputFileName)
		WriteVTIFile(self.Args.OutputFileName,ImageNew)

	
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will label the coronary CTA centerlines on the image using an extract centerline surface file")

        #Input filename of the coronary CTA
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="File name of the vti image containing the coronary CTA")

	#Input filename of the coronary segmented surface.
	parser.add_argument('-InputSurface', '--InputSurface', type=str, required=True, dest="InputSurface",help="File name of the input centerline fils containing the coronary centerlines")
	
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=True, dest="OutputFileName",help="File name in which to store the Image along with the Labels array.")
    
	args=parser.parse_args()
	ImageAnalysisLabelImage(args).Main()

