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

		print ("--- Create the Slices of the Myocardium")
		Slices,Centroid_slices=self.MyocardiumSlices(Volume,Centroid,Norm1,Size)

		print ("--- Morph the Volume to Match Surface Shape")
		self.MorphVolumeToSurface(Volume,Centroid,Norm1,Norm2,Apex,Size,Slices,Centroid_slices,Surface)	
	def MorphVolumeToSurface(self,Volume,Centroid,Norm1,Norm2,Apex,Size,Slices,Centroid_slices,Surface):
		#Get the number of points
		Npts=Volume.GetNumberOfPoints()	
		
		#Move all of the volume points to match surface shape
		for i in range(0,Npts):
			#Get the coordinate
			coord_=Volume.GetPoint(i)

			#Find the distance to myocardium centerline
			dist_P_to_line_=np.sqrt(vtk.vtkLine.DistanceToLine(coord_,Centroid,Apex))
			#Find the distance to myocardium apex
			dist_P_to_Apex_=np.power( np.power(coord_[0]-Apex[0],2) + np.power(coord_[1]-Apex[1],2) + np.power(coord_[2]-Apex[2],2),0.5)
                        #Compute the distance from Apex to Projected Point
			dist_Apex_to_ProjP_=np.power(np.power(dist_P_to_Apex_,2)-np.power(dist_P_to_line_,2),0.5)
                        ##Find the coordinate by moving from Apex a certain distance
			coord_ProjP_=Apex-Norm1*dist_Apex_to_ProjP_


			########### Create a line from LV Centerline to LV Mesh Point ########	
			#Get the Line Normals
			Norm_=np.sqrt(np.power(coord_[0]-coord_ProjP_[0],2)+np.power(coord_[1]-coord_ProjP_[1],2)+np.power(coord_[2]-coord_ProjP_[2],2))
			Norm_=(Norm_/np.linalg.norm(Norm_))
			#Create a vtkline
			line_=CreatePolyLine([coord_ProjP_,coord_ProjP_+Norm_*Size])
			
			########### Slice the Volume to Find Outer Most Point ##############
			centroid_,centroid_id=ClosestPoint(coord_,Centroid_slices)
			
			#Intersection between surface and line
			line_=CutPolyData(coord_,centroid_,Slices[centroid_id],Norm1)
			WriteVTPFile("Line_.vtp",line_)
			WriteVTPFile("Slice_.vtp",Slices[centroid_id])
			print (coord_)
			exit(1)	
        
	
	def MyocardiumSlices(self,Volume,Centroid,Norm1,Size):
                #Get the base and apex coordinate
		Coord_base=Centroid+Size*Norm1
		Coord_apex=Centroid-Size*Norm1
		#Create slice along the length of the appex
		Npts=int(((2*Size)/self.Args.SliceResolution))
		#Loop over the the myocardium and compute slices        
		Slices=[]
		Centroid_slices=[]
		progress_old=-1
		for i in range(Npts):
			progress_=self.PRINT_PROGRESS(i,Npts,progress_old)
			progress_old=progress_
			coord_=Coord_apex+Norm1*i*self.Args.SliceResolution
			slice_=ClippedSlices(coord_,Norm1,Volume)
			if slice_.GetNumberOfPoints()!=0:
				Slices.append(slice_)
				Centroid_slices.append(GetCentroid(slice_))
				WriteVTPFile("Slice%d.vtp"%i,slice_)
		return Slices,np.array(Centroid_slices)	

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
	
	parser.add_argument('-SliceResolution', '--SliceResolution', type=float, default=0.1, required=False, dest="SliceResolution",help="The number of slices for the myocardium.")

	args=parser.parse_args()
	ImageAnalysisMyocardiumMorphVolumeToSurface(args).Main()

		
				
			
			





	

