import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import numpy as np
import vtk
import argparse
from utilities import *
from PrincipleComponentAnalysis import PRINCIPLE_COMPONENT_ANALYSIS

class ImageAnalysisMyocardiumMorphVolumeToSurface():
	def __init__(self,Args):
		self.Args=Args

		if self.Args.OutputVolume is None:
			self.Args.OutputVolume=self.Args.InputVolume.replace(".vtu","_dilated.vtu")
	
	def Main(self):
                #Read the vtu file
		print ("--- Reading %s"%self.Args.InputVolume)
		Volume=ReadVTUFile(self.Args.InputVolume)
	
		#Read the vtp file of the LV from CTA
		print ("--- Read %s"%self.Args.InputSurface)
		Surface=ReadVTPFile(self.Args.InputSurface)

		#Case Rays from CTA Surface to Find Average of the CTA Volume
		print ("--- Project Volume Quantities onto Surface")
		self.ProjectVolumeToSurface(Volume,Surface)

	def ProjectVolumeToSurface(self,Volume,Surface):
		#Get PCA of the CTA Surface
		print ("--- Computing Principle Component Analysis of MBF Volume")
		Centroid,Norm1,Norm2,Apex,Size=PRINCIPLE_COMPONENT_ANALYSIS(Surface).main()
	
		#Get the number of points
		Npts=Surface.GetNumberOfPoints()
                
		#Distance from Apex To Centroid
		DistApexToCent=np.linalg.norm(Apex-Centroid)

		#Build a Cell Locator
		VolumeLocator = vtk.vtkCellLocator()
		VolumeLocator.SetDataSet(Volume)
		VolumeLocator.BuildLocator()
		
		#Move all of the volume points to match surface shape
		AverageMBF=np.zeros(Npts)
		Weights=np.zeros(Npts)
		for i in range(0,Npts):
			#Get the coordinate
			coord_=np.array(Surface.GetPoint(i))
			coord_ProjP_=ProjectedPointOnLine(coord_,Centroid,Apex,Norm1)
		
			#Get the Intersection of line with Volume Cells	
			cellids=vtk.vtkIdList()
			VolumeLocator.FindCellsAlongLine(coord_,coord_ProjP_,1e-6,cellids)
			
			#Get the Point Ids of the cells
			pointids=vtk.vtkIdList() #Point ids to store
			for j in range(cellids.GetNumberOfIds()):
				Volume.GetCellPoints(cellids.GetId(j), pointids)	 
			print(pointids.GetId(0))
			print(Volume.GetPoint(pointids.GetId(0)))
			print (cellids.GetId(0))
			exit(1)
				
	
			print (cellids.GetId(0))
			print(Volume.GetPoint(cellids.GetId(0)))
			exit(1)
			print (dir(Volume.GetPointData().GetArray("scalars")))
			exit(1)

			#Get the points from the cell ids
			for j in range(cellids.GetNumberOfIds()):
				print (cellids.GetId(j))
			exit(1)

			#Get the Distance from Apex to Point on Apex-Base Ais
			DistApexToCoordP_=np.linalg.norm(Apex-coord_ProjP_)

			#Project the Volumetric Points onto Surface. 
			#From Base to Centroid, Use Straight Line from Coord to Apex-Base Axis
			if DistApexToCoordP_>DistApexToCent:
				#Get the Norm
				norm_=np.array([coord_ProjP_[0]-coord_[0],coord_ProjP_[1]-coord_[1],coord_ProjP_[2]-coord_[2]])
				pSource=coord_ProjP_
			else: 
				norm_=np.array([Centroid[0]-coord_[0],Centroid[1]-coord_[1],Centroid[2]-coord_[2]])
				pSource=Centroid
				
			#Create a cut line through the dataset
			Slice_=CutPlane(Volume,coord_,Norm1)
			Line_ =CutLine(Slice_,coord_,coord_ProjP_,Norm1)
			
			#Get the Average Intensity
			MBF_=[Line_.GetPointData().GetArray("scalars").GetValue(j) for j in range(Line_.GetNumberOfPoints())]
			AverageMBF[i]=np.average(MBF_)
			Weights[i]=len(MBF_)
			WriteVTPFile("Slice_%d.vtp"%i,Slice_)	

		#Add Array to the Surface
		Surface=SurfaceAddArray(Surface,AverageMBF,"AverageMBF")
		Surface=SurfaceAddArray(Surface,Weights,"Weights")
	
		WriteVTPFile("Surface.vtp",Surface)
	
        #Print the progress of the loop
	def PRINT_PROGRESS(self,i,N,progress_old):
		progress_=(int((float(i)/N*100+0.5)))
		if progress_%10==0 and progress_%10!=progress_old: print ("    Progress: %d%%"%progress_)
		#progress_=int((float(i)/N)*1000)
		#if progress_%100==0: print ("    Progress: %d%%"%int(progress_/10.)),
		return progress_%10
	

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will dilate the myocardium by a specified vaue")
	parser.add_argument('-InputVolume', '--InputVolume', type=str, required=True, dest="InputVolume",help="The vtu file that contains the myocardial volume")
        
	parser.add_argument('-InputSurface', '--InputSurface', type=str, required=True, dest="InputSurface",help="The vtp file that contains the myocardium surface")

        #Output Filename 
	parser.add_argument('-OutputVolume', '--OutputVolume', type=str, required=False, dest="OutputVolume",help="The vtu file to store the output volume")
	
	args=parser.parse_args()
	ImageAnalysisMyocardiumMorphVolumeToSurface(args).Main()

		
				
			
			





	

