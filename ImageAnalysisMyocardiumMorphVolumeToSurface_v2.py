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
	
	def Main(self):
                #Read the vtu file
		print ("--- Reading %s"%self.Args.InputVolume)
		Volume=ReadVTUFile(self.Args.InputVolume)
	
		#Read the vtp file of the LV from CTA
		print ("--- Read %s"%self.Args.InputSurface)
		Surface=ReadVTPFile(self.Args.InputSurface)
		
		print ("--- Computing Princple Component Analysis and writing file")
		Centroid,Norm1,Norm2,Apex,Size=PRINCIPLE_COMPONENT_ANALYSIS(Volume).main()

		print ("--- Morph the Volume to Match Surface Shape")
		self.MorphVolumeToSurface(Volume,Centroid,Norm1,Norm2,Apex,Size,Surface)
	
	def MorphVolumeToSurface(self,Volume,Centroid,Norm1,Norm2,Apex,Size,Surface2):
		#Get the number of points
		Npts=Volume.GetNumberOfPoints()	

		#Get the outsurface of the MBF Volume
		Surface1=ExtractSurface(Volume)
		
		#Create OBB tree of Surfaces
		#Surface1=MBF Map
		#Surface2=Ventricle Surface from CTA
		obbTreeSurface1 = vtk.vtkOBBTree()
		obbTreeSurface1.SetDataSet(Surface1)
		obbTreeSurface1.BuildLocator()
		
		obbTreeSurface2 = vtk.vtkOBBTree()
		obbTreeSurface2.SetDataSet(Surface2)
		obbTreeSurface2.BuildLocator()

		pointsVTKintersection1 = vtk.vtkPoints()
		pointsVTKintersection2 = vtk.vtkPoints()

		#Move all of the volume points to match surface shape
		for i in range(0,Npts):
			#Get the coordinate
			coord_=np.array(Volume.GetPoint(i))
			norm_=np.array([coord_[0]-Centroid[0],coord_[1]-Centroid[1],coord_[2]-Centroid[2]])
			norm_=(norm_/np.linalg.norm(norm_))
			pTarget=coord_+norm_*Size/2.
			
			#Check intersection for VolumetricSurface
			code1= obbTreeSurface1.IntersectWithLine(Centroid,pTarget,pointsVTKintersection1, None)
			N1=pointsVTKintersection1.GetNumberOfPoints()
			
			#Check intersection for reference surface
			code2= obbTreeSurface2.IntersectWithLine(Centroid,pTarget,pointsVTKintersection2, None)
			N2=pointsVTKintersection2.GetNumberOfPoints()

			if N1>0 and N2>0:
				#Find the Furthest point from the centroid
				pointsIntersection1=np.array([pointsVTKintersection1.GetPoint(j) for j in range(N1)])
				coord1_,coord1id_=FurthestPoint(Centroid,pointsIntersection1)
	
				#Get the coordinates of Volume surface and CTA surface
				coord2_=pointsVTKintersection1.GetPoint(N2-1)

				#Move the coordinate in normal direction
				norm1_=np.array([coord1_[0]-coord2_[0],coord1_[1]-coord2_[1],coord1_[2]-coord2_[2]])
				Volume.GetPoints().SetPoint(i,coord_-norm1_)
			else:
				Volume.GetPoints().SetPoint(i,coord_)	
		WriteVTUFile("abc2.vtu",Volume)
		exit(1)

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
	parser.add_argument('-OutputVolume', '--OutputVolume', type=str, required=True, dest="OutputVolume",help="The vtu file to store the output volume")
	
	args=parser.parse_args()
	ImageAnalysisMyocardiumMorphVolumeToSurface(args).Main()

		
				
			
			





	

