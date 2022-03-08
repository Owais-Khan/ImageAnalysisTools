import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *
from PrincipleComponentAnalysis import PRINCIPLE_COMPONENT_ANALYSIS

class ImageAnalysisMyocardiumMorphVolumeToSurface():
	def __init__(self,Args):
		self.Args=Args

		if self.Args.OutputSurface is None:
			self.Args.OutputSurface=self.Args.InputSurface.replace(".vtp","_ProjectedMBF.vtp")
			self.Args.OutputThickness=self.Args.InputSurface.replace(".vtp","_ProjectedMBFThicknessMap.dat")
		else:
			self.Args.OutputThickness=self.Args.OutputSurface.replace(".vtp","_ThicknessMap.dat")
			
	def Main(self):
                #Read the vtu file
		print ("--- Reading %s"%self.Args.InputVolume)
		
		if   self.Args.InputVolume[-4:]==".vtk":
			Volume=ReadVTKFile(self.Args.InputVolume)   #Read a VTK volume stack
			Volume=ThresholdByUpper(Volume,"scalars",0) #Convert to an unstructured grid 
		
		elif self.Args.InputVolume[-4:]==".vtu":
			Volume=ReadVTUFile(self.Args.InputVolume) #Read a VTU unstructured grid
		else:
			print ("The extension %s is not valid for volume"%self.Args.InputVolume[-4:])
			print ("Exiting...")
			exit(1)

	
		#Read the vtp file of the LV from CTA
		print ("--- Read %s"%self.Args.InputSurface)
		Surface=ReadVTPFile(self.Args.InputSurface)

		#Case Rays from CTA Surface to Find Average of the CTA Volume
		print ("--- Project Volume Quantities onto Surface")
		self.ProjectVolumeToSurface(Volume,Surface)

	def ProjectVolumeToSurface(self,Volume,Surface):
		#Get PCA of the CTA Surface
		print ("------ Computing Principle Component Analysis of CTA Surface")
		Centroid,Norm1,Norm2,Apex,Size=PRINCIPLE_COMPONENT_ANALYSIS(Surface).main()
	
		print ("------ Computing Surface Normals of CTA Surface")
		Surface=SurfaceNormals(Surface)

		#Get the number of points
		Npts=Surface.GetNumberOfPoints()
                
		#Build a Cell Locator
		VolumeLocator = vtk.vtkCellLocator()
		VolumeLocator.SetDataSet(Volume)
		VolumeLocator.BuildLocator()
	
		#Create a dictionary to store Average, Median and 75th Percentile MBF
		MBF={"MBF_WallAveraged":np.zeros(Npts), "Weights":np.zeros(Npts)}
	
		#Create an array to store surface points that have zero values
		NodesZeroMBF=[]

		FileOutputThickness=open(self.Args.OutputThickness,'w')
		FileOutputThickness.write("NormalizedDistanceFromApex Npts DistanceFromOutWall MBF \n")
		for i in range(0,Npts):
			coord_=np.array(Surface.GetPoint(i))

			#Get the Coordinate Projected on Base Apex Line
			coord_ProjP_=ProjectedPointOnLine(coord_,Centroid,Apex,Norm1)
			#Distance from Apex
			DistToApex_=np.linalg.norm(coord_ProjP_-Apex)

			#Get the Surface Normal
			CellNormal_=np.array(Surface.GetPointData().GetArray("Normals").GetTuple(i))
			
			#Get Source and Target Points in Normal Directions
			DistToCenter=np.linalg.norm(coord_-Centroid)
			pSource=coord_-CellNormal_*DistToCenter*0.75
			pTarget=coord_+CellNormal_*DistToCenter*0.75

			#Get the Intersection of line with Volume Cells	
			cellids=vtk.vtkIdList()
			VolumeLocator.FindCellsAlongLine(pSource,pTarget,1e-6,cellids)
			
			#Get the Point Ids of the cells that intersect with line
			pointids=vtk.vtkIdList() #Point ids to store
			for j in range(cellids.GetNumberOfIds()):
				pointids_=vtk.vtkIdList()
				Volume.GetCellPoints(cellids.GetId(j), pointids_)
				#Add all four points of the hex cell into the point id list
				for k in range(pointids_.GetNumberOfIds()):
					pointids.InsertNextId(pointids_.GetId(k))

			#Get the MBF Data using the point ids
			MBF_=np.array([Volume.GetPointData().GetArray("scalars").GetValue(pointids.GetId(j)) for j in range(pointids.GetNumberOfIds())])

			#Get the distance from the outer wall to calculate thickness
			DistFromOutWallArray_=np.array([np.linalg.norm(pTarget-Volume.GetPoint(pointids.GetId(j))) for j in range(pointids.GetNumberOfIds())])
			DistSortedIds_=np.argsort(DistFromOutWallArray_)
	
			if len(MBF_)>0:
				#Get the average MBF values of all the point nodes
				MBF["MBF_WallAveraged"][i]=np.average(MBF_)
				MBF["Weights"][i]=len(MBF_)
				
				#Write the thickness values to output file
				#Calculate Thickness from the Outer wall
				FileOutputThickness.write("%.05f %d "%(DistToApex_/Size,len(MBF_)))
				FileOutputThickness.write("0.0 %.05f "%MBF_[DistSortedIds_[0]])
				for j in range(1,len(MBF_)):
					DistFromOutWall_=DistFromOutWallArray_[DistSortedIds_[j]]
					DistFromOutWall_-=DistFromOutWallArray_[DistSortedIds_[0]]
					FileOutputThickness.write("%.05f %.05f "%(DistFromOutWall_,MBF_[DistSortedIds_[j]]))
				FileOutputThickness.write("\n")	

			else:
                                
				MBF["MBF_WallAveraged"][i]=-999.0
				MBF["Weights"][i]=-999.0 
				NodesZeroMBF.append(i)				
		
		#Close the File
		FileOutputThickness.close()	

		#Fill all of the empty values with values from average of the neightbouring cells.
		print ("------ There are %d points that have Zero MBF value"%len(NodesZeroMBF))
		print ("------ Filling these values with MBF from Neighbouring points")
		for PointId in NodesZeroMBF:
			coord_=Surface.GetPoint(PointId)
			Dist_=np.array([np.linalg.norm(coord_-np.array(Surface.GetPoint(i))) for i in range(Npts)])
			DistSortIds_=np.argsort(Dist_)
			MBF_closest_=[MBF["MBF_WallAveraged"][i] for i in DistSortIds_]
			MBF_closest_=[Value for Value in MBF_closest_ if Value!=-999.0]
			MBF["MBF_WallAveraged"][PointId]=np.average(MBF_closest_[0:5])
			MBF["Weights"][PointId]=1
			

		print ("--- Computing 50th and 75th percentiles")	
		#Get 75th Percentile MBF
		MBF_75Q=np.percentile(MBF["MBF_WallAveraged"],0.75)
		MBF_50Q=np.percentile(MBF["MBF_WallAveraged"],0.50)

		#Add Array to the Surface
		Surface=SurfaceAddArray(Surface,MBF["MBF_WallAveraged"],"MBF_WallAveraged")
		Surface=SurfaceAddArray(Surface,MBF["MBF_WallAveraged"]/MBF_75Q,"MBF_Normalized75Q")
		Surface=SurfaceAddArray(Surface,MBF["MBF_WallAveraged"]/MBF_50Q,"MBF_Normalized50Q")
		Surface=SurfaceAddArray(Surface,MBF["Weights"],"Weights")

		WriteVTPFile(self.Args.OutputSurface,Surface)
	
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
	parser.add_argument('-OutputSurface', '--OutputSurface', type=str, required=False, dest="OutputSurface",help="The vtp file to store the output surface with MBF projected from volume")
	
	args=parser.parse_args()
	ImageAnalysisMyocardiumMorphVolumeToSurface(args).Main()

		
				
			
			





	

