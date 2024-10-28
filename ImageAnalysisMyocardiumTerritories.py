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

		if self.Args.OutputFileName is None:
			if self.Args.InputFileName.find("/")>=0: self.Args.OutputFilename=self.Args.InputFileName.replace(self.Args.InputFileName.split("/")[-1],"MBF_Territories.vtu")
			else: self.Args.OutputFileName="MBF_Territories.vtu"
			self.Args.OutputFilename2=self.Args.OutputFileName.replace(".vtu","_Labels.dat")

		#Coronary Centerline Files
		LADFiles        =sorted(glob("%s/mesh-surfaces/wall_L_LAD_*.vtp"%self.Args.VesselSurfaces))
		IntermediusFiles=sorted(glob("%s/mesh-surfaces/wall_L_Intermedius_*.vtp"%self.Args.VesselSurfaces))
		LCxFiles        =sorted(glob("%s/mesh-surfaces/wall_L_LCx_*.vtp"%self.Args.VesselSurfaces))
		Diag1Files      =sorted(glob("%s/mesh-surfaces/wall_L_Diag1_*.vtp"%self.Args.VesselSurfaces))
		Diag2Files      =sorted(glob("%s/mesh-surfaces/wall_L_Diag2_*.vtp"%self.Args.VesselSurfaces))
		PDAFiles        =sorted(glob("%s/mesh-surfaces/wall_R_PDA_*.vtp"%self.Args.VesselSurfaces))
		PLAFiles        =sorted(glob("%s/mesh-surfaces/wall_R_PLA_*.vtp"%self.Args.VesselSurfaces))
		
		#Merge all of the files
		self.CenterlineFileNamesPaths=LADFiles+IntermediusFiles+LCxFiles+Diag1Files+Diag2Files+PDAFiles+PLAFiles
	
		self.CenterlineFileNames=[filename.split("/")[-1] for filename in self.CenterlineFileNamesPaths] 

		count=0
		for CenterlineFileName_ in self.CenterlineFileNames:
			print ("Surface %d is: %s"%(count,CenterlineFileName_))
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


		#Find the array name with work "scalars"
		if self.Args.ArrayName is None:
			for i in range(VentricleMesh.GetPointData().GetNumberOfArrays()):
				ArrayName_=str(VentricleMesh.GetPointData().GetArrayName(i))
				if ArrayName_.find("calars")>0: self.Args.ArrayName=ArrayName_
		print ("--- The array that contains MBF values is: %s"%self.Args.ArrayName)

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
		print ("--- Writing Territory Labels to %s"%self.Args.OutputFilename2)
		outfile=open(self.Args.OutputFilename2,'w')
		outfile.write("TerritoryLabel CenterlineName\n")
		counter=0
		for filename in filenames:
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
			outfile.write("%d %s\n"%(counter,filename.split("/")[-1]))
			counter+=1 
		outfile.close()
		return CL_coords	

	def Get_Point_Data(self,LV_mesh):
		print ("--- Getting LV Coordinates and MBF")
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
	parser = argparse.ArgumentParser(description="This script will compute terriorites of the myocardium using coronary vessel centerlines/surfaces.")

        #Input filename of the perfusion map
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="Volumetric MBF data in VTK or VTU format")
	
	#Input folder that contains the coronary centerlines
	parser.add_argument('-VesselSurfaces', '--VesselSurfaces', type=str, required=True, dest="VesselSurfaces",help="Folder that contains the centerline files labels into L_*.vtp and R_*.vtp tags.")
	
        #Array Name of the Data
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName",help="The name of the array containing the MBF values.")

        #Output argumenets
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename of the volumetric data with territories")
        
	args=parser.parse_args()
	ImageAnalysisMyocardiumCoronaryTerritories(args).main()

