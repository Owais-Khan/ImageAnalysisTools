import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
import argparse
from utilities import *

class ImageAnalysisMyocardiumCoronaryTerritories():
	def __init__(self,Args):
		#Input argumenets
		self.Args=Args

		#Coronary Centerline Files
		self.CenterlineFileNamesPaths=sorted(glob("%s/*.vtp"%self.Args.CenterlinesFolder))
		self.CenterlineFileNames=[filename.split("/")[-1] for filename in self.CenterlineFileNamesPaths] 
              
		#Separate the Coronary Territories into Left and Right side 
		"""self.CenterlineFileNamesPaths=[]
		self.CenterlineFileNames=[]
		for i in range(len(CenterlineFileNames)):
			if CenterlineFileNames[i].find("L_")>=0:
				self.CenterlineFileNames.append(CenterlineFileNames[i])
				self.CenterlineFileNamesPaths.append(CenterlineFileNamesPaths[i])
			if CenterlineFileNames[i].find("R_")>=0:
				self.CenterlineFileNames.append(CenterlineFileNames[i])
				self.CenterlineFileNamesPaths.append(CenterlineFileNamesPaths[i])"""

		count=0
		for CenterlineFileName_ in self.CenterlineFileNames:
			print ("Centerline %d is: %s"%(count,CenterlineFileName_))
			count+=1

	def main(self):
		#Get the areas for all of the surface caps
		#Centroids=self.Get_Centroids(self.CenterlineFileNamesPaths,self.CenterlineFileNames)

                #Read the vtu file
		print ("--- Reading Left Ventricle Myocardial Blood Flow Data: %s"%self.Args.InputFileName)
		if   self.Args.InputFileName[-4:]==".vtk":
			VentricleMesh=ReadVTKFile(self.Args.InputFileName)   #Read a VTK volume stack
			VentricleMesh=ThresholdByUpper(VentricleMesh,self.Args.ArrayName,1) #unstruct grid 
		elif self.Args.InputFileName[-4:]==".vtu":
			VentricleMesh=ReadVTUFile(self.Args.InputFileName) #Read a VTU unstruct grid
		else:
			print ("The extension %s is not valid for volume"%self.Args.InputFileName[-4:])
			print ("Exiting...")
			exit(1)



		#Get the Nodes and Scalars of the LV Mesh
		print ("--- Getting Coordinates from the Left Ventricle Mesh")
		VentricleCoordinates,ImageScalars=self.Get_Point_Data(VentricleMesh)

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
		writer=vtk.vtkXMLUnstructuredGridWriter()
		writer.SetInputData(LV_mesh)
		if self.Args.OutputFileName is None:
			self.Args.OutputFilename=self.InputFileName.replace(self.InputputFileName.split("/")[-1],"ImageAnalysisMyocardiumTerritories.vtu")
		writer.SetFileName(self.Args.OutputFileName)
		writer.Update()
			
	
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

	def Get_Point_Data(self,LV_mesh):
		print ("------ Getting LV Coordinates and MBF")
		N=LV_mesh.GetNumberOfPoints()	
		Coords=np.zeros(shape=(N,3))
		imageScalars=np.zeros(N)
		for i in range(N):
			coord_=LV_mesh.GetPoint(i)
			Coords[i,0]=coord_[0]
			Coords[i,1]=coord_[1]
			Coords[i,2]=coord_[2]
			imageScalars[i]=LV_mesh.GetPointData().GetArray(self.Args.ArrayName).GetValue(i)
		return Coords,imageScalars	

	def Read_Vtu(self,filename):
		reader=vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(filename)
		reader.Update()
		return reader.GetOutput()				
		

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
	parser = argparse.ArgumentParser(description="This script will compute terriorites of the myocardium using coronary vessel centerlines.")

        #Input filename of the perfusion map
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="Volumetric MBF data in VTK or VTU format")
	
	#Input folder that contains the coronary centerlines
	parser.add_argument('-CenterlinesFolder', '--CenterlinesFolder', type=str, required=True, dest="CenterlinesFolder",help="Folder that contains the centerline files labels into L_*.vtp and R_*.vtp tags.")
	
        #Array Name of the Data
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False,default="scalars", dest="ArrayName",help="The name of the array containing the MBF values.")

        #Output argumenets
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename of the volumetric data with territories")
        
	args=parser.parse_args()
	ImageAnalysisMyocardiumCoronaryTerritories(args).main()

