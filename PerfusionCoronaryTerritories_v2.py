import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
import argparse
from utilities import *

class ImageAnalysisPerfusionCoronaryTerritories():
	def __init__(self,Args):
		#Input argumenets
		self.Args=Args
		
		if self.Args.InputFileName is None:
			self.Args.InputFileName="Results/MBF_ProjectedToSurface.vtp"

		if self.Args.CenterlinesFolder is None:
			self.Args.CenterlinesFolder="Results/Centerlines/"

		if self.Args.OutputFileName is None:
                        self.Args.OutputFileName="Results/MBF_ProjectedToSurface_Territories.vtp"

		#The LV filename that contains MBF
		self.VentricleFileName=self.Args.InputFileName

		#Coronary Centerline Files
		CenterlineFileNamesPaths=sorted(glob("%s/wall_*cl.vtp"%self.Args.CenterlinesFolder))
		CenterlineFileNames=[filename.split("/")[-1] for filename in CenterlineFileNamesPaths] 
              
		#Separate the Coronary Territories into Left and Right side 
		self.CenterlineFileNamesPaths=[]
		self.CenterlineFileNames=[]
		for i in range(len(CenterlineFileNames)):
			CLName_=CenterlineFileNames[i]
			if CLName_.find("RCA")>=0 or CLName_.find("LCx")>=0 or CLName_.find("LAD")>=0 or CLName_.find("Ramus")>=0 or CLName_.find("Diag")>=0:
				self.CenterlineFileNames.append(CenterlineFileNames[i])
				self.CenterlineFileNamesPaths.append(CenterlineFileNamesPaths[i])


	def Main(self):
		#Get the areas for all of the surface caps
		Centroids=self.Get_Centroids(self.CenterlineFileNamesPaths,self.CenterlineFileNames)

		#Read the LV mesh file
		print ("--- Reading Left Ventricle Myocardial Blood Flow Data: %s"%self.VentricleFileName)
		VentricleMesh=ReadVTPFile(self.VentricleFileName)

		#Get the Nodes and Scalars of the LV Mesh
		print ("--- Getting Coordinates from the Left Ventricle Mesh")
		VentricleCoordinates=self.GetPointData(VentricleMesh)

		#Read the Centerline Coordinates
		print ("--- Reading Centerline Files")
		CenterlineCoordinates=self.Read_CL_Coords(self.CenterlineFileNamesPaths)
		
		#For each LV node, find the shortest distance to the CL
		print ("--- Getting the shortest distance for each Ventricle Node to closest centerline")
		ClosestCenterLines=self.Get_Shortest_Distance(CenterlineCoordinates,VentricleCoordinates)	
		#Now write the data with a new territory division
		print ("--- Writing the territory maps")
		self.Write_Territories(VentricleMesh,ClosestCenterLines)
	
	def Write_Territories(self,LV_mesh,Closest_CL):
		#Create New Array in LV mesh with CL names
		territories=vtk.vtkIntArray()
		territories.SetNumberOfComponents(1)
		territories.SetName("TerritoryMaps")
		N=LV_mesh.GetNumberOfPoints()
		for i in range(N):
			index=self.CenterlineFileNames.index(Closest_CL[i])
			territories.InsertNextValue(index)
                
		LV_mesh.GetPointData().AddArray(territories)
		WriteVTPFile(self.Args.OutputFileName,LV_mesh)
	
	def Get_Shortest_Distance(self,CL_coords,LV_coords):
		#For each Node, Loop over all of the centerlines
		N=len(LV_coords)
		Closest_CL=[]
		count_p=0.05
		progress_old=-1
		for i,LV_coord in enumerate(LV_coords):
			progress_=self.PRINT_PROGRESS(i,N,progress_old)
			progress_old=progress_
			distance={}
			#Loop over all of the centerlines
			for Key in CL_coords.keys():
				distance[Key]=np.min(DISTANCE.cdist(CL_coords[Key],[LV_coords[i]]))
			#Store the Centerline that is closest to the point
			Closest_CL.append(min(distance,key=distance.get))
		return Closest_CL	

	def Read_CL_Coords(self,filenames):
		CL_coords={}
		for filename in filenames:
			print ("-----%s"%filename)
			filename_short=filename.split("/")[-1]
			reader=vtk.vtkXMLPolyDataReader()
			reader.SetFileName(filename)
			reader.Update()
			CL_data=reader.GetOutput()
			N=CL_data.GetNumberOfPoints()
			points_=np.zeros(shape=(N,3))
			for i in range(N):
				points_[i]=CL_data.GetPoint(i)	
			CL_coords[filename_short]=points_
		return CL_coords	

	def GetPointData(self,LV_mesh):
		print ("------ Getting LV Coordinates and MBF")
		N=LV_mesh.GetNumberOfPoints()	
		Coords=np.zeros(shape=(N,3))
		imageScalars=np.zeros(N)
		for i in range(N):
			coord_=LV_mesh.GetPoint(i)
			Coords[i,0]=coord_[0]
			Coords[i,1]=coord_[1]
			Coords[i,2]=coord_[2]
		return Coords	

	def Get_Centroids(self,surface_files,filenames):
		print ("-"*50)
		print ("------ Computing Centroids of Surface Outlets")
		
		Centroids={}
		count=0
		for surface_file in surface_files:
			print ("----Processing %s"%filenames[count])
			#Read the Surface File
			reader=vtk.vtkXMLPolyDataReader()
			reader.SetFileName(surface_file)
			reader.Update()
			
			#Get the Centroid
			centerOfMassFilter=vtk.vtkCenterOfMass()
			centerOfMassFilter.SetInputData(reader.GetOutput())
			centerOfMassFilter.SetUseScalarsAsWeights(False)
			centerOfMassFilter.Update()
			center_=centerOfMassFilter.GetCenter()
			Centroids[filenames[count]]=[center_[0],center_[1],center_[2]]
			del reader,centerOfMassFilter
			count+=1
		return Centroids

        
	def PRINT_PROGRESS(self,i,N,progress_old):
		progress_=(int((float(i)/N*100+0.5)))
		if progress_%10==0 and progress_%10!=progress_old: print ("    Progress: %d%%"%progress_)
		return progress_%10



if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute terriorites of the myocardium using coronary vessels.")

        #Input filename of the perfusion map
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=False, dest="InputFileName",help="Volumetric Mesh that contains the Myocardial Blood Flow Data")
	
	#Input folder that contains the coronary centerlines
	parser.add_argument('-CenterlinesFolder', '--CenterlinesFolder', type=str, required=False, dest="CenterlinesFolder",help="Folder that contains the centerline files")
	
        #Array Name of the Data
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False,default="ImageScalars", dest="ArrayName",help="The name of the array containing the MBF values")

        #Output argumenets
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename of the volumetric data with territories")
        
	args=parser.parse_args()
	ImageAnalysisPerfusionCoronaryTerritories(args).Main()

