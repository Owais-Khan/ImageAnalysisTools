import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import os
from glob import glob
import argparse
from utilities import *
import matplotlib.pyplot as plt
class ImageAnalysisProjectImageWallToSurface():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFileName is None:
			#self.Args.OutputFileName=self.Args.InputSurface[0:-4]+"_"+self.Args.InputImage.split("/")[-1][0:-4]+self.Args.InputSurface[-4:]
			self.Args.OutputFileName=self.Args.InputImage.replace(".vti","_Surface.vtp")
			print ("--- OutputFileName not provided...")
			print ("--- Using %s"%self.Args.OutputFileName)
			
	def Main(self):
		#Read the source data
		print ("--- Loading the source data: %s"%self.Args.InputImage)
		if self.Args.InputImage[-4:]==".vti":
			SourceData=ReadVTIFile(self.Args.InputImage)
		else:
			print ("Input file format not detected. Exiting...")
			exit(1)

		print ("--- Loading the target volume/surface mesh: %s"%self.Args.InputSurface)
		if self.Args.InputSurface[-4:]==".vtp":
			TargetData=ReadVTPFile(self.Args.InputSurface)
		else:
			print ("Target surface format not detected. Exiting...")
			exit(1)	

		#Convert Image Point Data to CellData
		ConvertToCellData=vtk.vtkPointDataToCellData()
		ConvertToCellData.SetInputData(SourceData)
		ConvertToCellData.Update()
		SourceData=ConvertToCellData.GetOutput()		


		#Compute the Normal Vectors
		print ("Computing Surface Normals ...")
		ComputeNormals_=vtk.vtkPolyDataNormals()
		ComputeNormals_.SetInputData(TargetData)
		ComputeNormals_.ComputePointNormalsOn()
		ComputeNormals_.ComputeCellNormalsOff()
		ComputeNormals_.Update()
		TargetData=ComputeNormals_.GetOutput()


                #Build a Cell Locator
		print ("Building a volume locator...")
		VolumeLocator = vtk.vtkCellLocator()
		VolumeLocator.SetDataSet(SourceData)
		VolumeLocator.BuildLocator()

		#Create an array for signal intensities
		SignalIntensitiesAvg=np.zeros(TargetData.GetNumberOfPoints())
		SignalIntensitiesMax=np.zeros(TargetData.GetNumberOfPoints())

		#Now Project the Image Intensities onto the Surface
		print ("Looping over the surface points ...")
		for i in range(0,TargetData.GetNumberOfPoints()):
			#Get the normal vectors
			Normals_=np.array(TargetData.GetPointData().GetArray("Normals").GetTuple(i))
                        
			#Get Source and Target Points in Normal Directions
			pSource=TargetData.GetPoint(i)+Normals_*self.Args.Thickness #Pointing outwards
			pTarget=TargetData.GetPoint(i)-Normals_*self.Args.Thickness #Pointing inwards

			#Create a line
			Line=ConvertPointsToLine(np.array([pSource,pTarget]))
			
                        #Get the Intersection of line with Volume Cells 
			IntersectCellIds_=vtk.vtkIdList()
			VolumeLocator.FindCellsAlongLine(pSource,pTarget,1e-6,IntersectCellIds_)

			SignalIntensities_=np.zeros(IntersectCellIds_.GetNumberOfIds())
			for j in range(IntersectCellIds_.GetNumberOfIds()):
				id_=IntersectCellIds_.GetId(j)
				value_=SourceData.GetCellData().GetArray(self.Args.ArrayName).GetValue(id_)
				SignalIntensities_[j]=value_
				
			SignalIntensitiesAvg[i]=np.average(SignalIntensities_)
			SignalIntensitiesMax[i]=np.max(SignalIntensities_)

		print ("The Average Signal Intensity is: %.05f"%np.average(SignalIntensitiesAvg))
		print ("The Median  Signal Intensity is: %.05f"%np.percentile(SignalIntensitiesAvg,50))
		print ("The standard deviation is        %.05f"%np.std(SignalIntensitiesAvg))
		TargetData=SurfaceAddArray(TargetData,SignalIntensitiesAvg,"SignalIntensityAvg")
		TargetData=SurfaceAddArray(TargetData,SignalIntensitiesMax,"SignalIntensityMax")

		print ("Writing OutputFile: %s"%self.Args.OutputFileName)
		WriteVTPFile(self.Args.OutputFileName,TargetData)	

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will extract the signal intensities from the image wall and project it onto the surface for visualization. There is option to choose two methods: i) Full-width half max (FWHM) or ii) constant thickness. It is assumed that image and surface mesh are registered adequately.")

	parser.add_argument('-InputImage', '--InputImage', type=str, required=True, dest="InputImage",help="File name of the image that contains the signal intensities.")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName",default="ImageScalars",help="The array name that contains the image intensitites.")
        
	#Input filename of the coronary segmented surface.
	parser.add_argument('-InputSurface', '--InputSurface', type=str, required=True, dest="InputSurface",help="File name of the input surface onto which the signal intensities are to be projcted. ")
	
	parser.add_argument('-Thickness', '--Thickness', type=float, required=False, default=1.0, dest="Thickness",help="File thickness of the wall over which to average the signal intensities. ")
        
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="File name in which to store the output data.")
        
	args=parser.parse_args()
       
	ImageAnalysisProjectImageWallToSurface(args).Main()

