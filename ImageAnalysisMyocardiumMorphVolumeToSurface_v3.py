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
		
		print ("--- Computing Princple Component Analysis and writing file")
		Centroid,Norm1,Norm2,Apex,Size=PRINCIPLE_COMPONENT_ANALYSIS(Volume).main()
		
		print ("--- Computing Principle Component Analysis of LV Surface")	

		print ("--- Morph the Volume to Match Surface Shape")
		self.MorphVolumeToSurface(Volume,Centroid,Norm1,Norm2,Apex,Size,Surface)
	
	def MorphVolumeToSurface(self,Volume,Centroid,Norm1,Norm2,Apex,Size,Surface):
		#Compute the PCA of LV Surface
		CentroidSurf,Norm1Surf,Norm2Surf,ApexSurf,SizeSurf=PRINCIPLE_COMPONENT_ANALYSIS(Surface).main()

		#Get the number of points
		NptsV=Volume.GetNumberOfPoints()
		VolumeSurface=ExtractSurface(Volume)

                #Create OBB tree of Surfaces
		obbTreeSurface1 = vtk.vtkOBBTree()
		obbTreeSurface1.SetDataSet(VolumeSurface)
		obbTreeSurface1.BuildLocator()
		pointsVTKintersection = vtk.vtkPoints()
	
		NptsS=Surface.GetNumberOfPoints()
		SurfacePoints=np.array([Surface.GetPoint(i) for i in range(NptsS)])
		
		#Move all of the volume points to match surface shape
		for i in range(0,NptsV):
			#Get the coordinate
			coord_=np.array(Volume.GetPoint(i))

			#Find the location (coord,distance) on the LV Apex-Base axis
			dist_P_to_line_=np.sqrt(vtk.vtkLine.DistanceToLine(coord_,Centroid,Apex))
			dist_P_to_Apex_=np.power( np.power(coord_[0]-Apex[0],2) + np.power(coord_[1]-Apex[1],2) + np.power(coord_[2]-Apex[2],2),0.5)
			dist_Apex_to_ProjP_=np.power(np.power(dist_P_to_Apex_,2)-np.power(dist_P_to_line_,2),0.5)
			coord_ProjP_=Apex+Norm1*dist_Apex_to_ProjP_

			#Get the Norm
			norm_=np.array([Centroid[0]-coord_[0],Centroid[1]-coord_[1],Centroid[2]-coord_[2]])
			norm_=norm_/np.linalg.norm(norm_)
			
			#Find the closest coordinate on LV Surface
			coord2_,coord2_id_=ClosestPoint(coord_,SurfacePoints)		
			
			#Find the distance from Centroid to MBF and Surface Coordinates
			distCtoV_=np.linalg.norm(coord_-Centroid)
			distCtoS_=np.linalg.norm(coord2_-Centroid)
		
			#Check intersection for VolumetricSurface
			pTarget=coord_-norm_*Size #Volume 
			code=obbTreeSurface1.IntersectWithLine(Centroid,pTarget,pointsVTKintersection, None)
			N=pointsVTKintersection.GetNumberOfPoints()
			pointsIntersection=np.array([pointsVTKintersection.GetPoint(j) for j in range(N)])

			#if len(pointsIntersection)==0: coord1_=coord_ #If no intersection, keep same coord
			print (pointsIntersection)
			print (coord_)
			coord1_,coord1id_=FurthestPoint(Centroid,pointsIntersection)

			#Find the distance between surface and volume
			dist_=np.linalg.norm(coord1_-coord2_)	
			
			#Strech the coordinates
			if distCtoV_<=distCtoS_: coord_new_=coord_-norm_*dist_#*(1-dist_Apex_to_ProjP_/Size)
			else: coord_new_=coord_+norm_*dist_#*(1-dist_Apex_to_ProjP_/Size)
	
			#Update the volume
			Volume.GetPoints().SetPoint(i,coord_new_)

		WriteVTUFile("abc2.vtu",Volume)

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

		
				
			
			





	

