#This script will read the volumetric mesh
#and remove all scalar values close to the
#surface boundary by the given number of layers


import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
import argparse

class ImageAnalysisMyocardiumRemoveBoundaryData():
	def __init__(self,Args):
		self.Args=Args

		#If output filename is not defined, create your own
		if self.Args.OutputFileName is None:
			InputFileNameStripped=self.Args.InputFileName.split("/")[-1]
			OutputFolder =self.Args.InputFileName.replace(InputFileNameStripped,"")
			OutputFileName=OutputFolder+"MyocardiumRemoveBoundaryData"
			self.Args.OutputFileName=OutputFileName

	def main(self):
		#Read the LV mesh file that contains the territory mappings and MBF
		print ("--- Reading Left Ventricle Myocardial Blood Flow Data and Territory Maps: %s"%self.Args.InputFileName)
		Volume=self.Read_Vtu(self.Args.InputFileName)

		#Now extract all of the surface mesh nodes fromt the volumetric data
		print ("--- Extracting the Surface Coordinates from the Volumetric Mesh")
		SurfaceCoords=self.ExtractSurfaceCoordinates(Volume)

		#Now Find the distance from the each of the volumetric node to the surface coordinate
		print ("--- Scrap off N=%d layers from the Volume"%self.Args.Layers)
		ThresholdArray=[False]*Volume.GetNumberOfPoints()
		ThresholdArray=self.ProximityToSurface(Volume,SurfaceCoords,ThresholdArray,0)
		
	def ThresholdVolume(self,Volume,ThresholdArray):
		N=Volume.GetNumberOfPoints()
		for i in range(N):
			if ThresholdArray[i] is True:
				Volume.GetPointData().GetArray(self.Args.ArrayName).SetValue(i,0)
		return Volume	

	def ProximityToSurface(self,Volume,SurfaceCoords,ThresholdArray,Counter):
		if Counter==self.Args.Layers:
			print ("Finished Removing %d Layers from the Myocardium walls"%self.Args.Layers)
			return ThresholdArray
 
		Np_vol =Volume.GetNumberOfPoints()
		Np_surf=len(SurfaceCoords[0])
		progress_old=-1	
		
		#Create a list of Surface Nodes to first check overlap with Volume node
		#Then check the distance. This may be faster 
		SurfaceCoordsX_aslist=np.ndarray.tolist(SurfaceCoords[0])
		SurfaceCoordsY_aslist=np.ndarray.tolist(SurfaceCoords[1])
		SurfaceCoordsZ_aslist=np.ndarray.tolist(SurfaceCoords[2])
		SurfaceCoords_aslist=[SurfaceCoordsX_aslist,SurfaceCoordsY_aslist,SurfaceCoordsZ_aslist]

		#Loop over all of the volumetric nodes
		for i in range(Np_vol):
			#print the progress so far
			progress_=self.PRINT_PROGRESS(i,Np_vol,progress_old)
			progress_old=progress_
		
			#Loop over the volume points and shed the layer	
			volume_coord_=Volume.GetPoint(i)
			if (ThresholdArray[i] is False) and (volume_coord_[0] in SurfaceCoords_aslist[0]) and (volume_coord_[1] in SurfaceCoords_aslist[1]) and (volume_coord_[2] in SurfaceCoords_aslist[2]):
				ThresholdArray[i]=True
				SurfaceCoords_aslist[0].remove(volume_coord_[0])
				SurfaceCoords_aslist[1].remove(volume_coord_[1])
				SurfaceCoords_aslist[2].remove(volume_coord_[2])

                
		#Now Make all of the scalar values close to the surface 0.
		Volume=self.ThresholdVolume(Volume,ThresholdArray)
		self.WriteOutput(self.Args.OutputFileName,Volume)	
		
		#Get New surface Coordinates
		Threshold=vtk.vtkThreshold()
		Threshold.SetInputData(Volume)
		Threshold.ThresholdByUpper(1)
		Threshold.SetInputArrayToProcess(0, 0, 0,"vtkDataObject::FIELD_ASSOCIATION_POINTS",self.Args.ArrayName)
		Threshold.Update()
		Threshold=Threshold.GetOutput()
		SurfaceCoords=self.ExtractSurfaceCoordinates(Threshold)
		
		#Save the current volume with "i"th layer scrapped off
		print ("--- Writing File: %s"%self.Args.OutputFileName+"_layer_%d.vtu"%Counter)
		self.WriteOutput(self.Args.OutputFileName+"_layer_%d.vtu"%Counter,Volume)
	
		#Repeat this in a loop until layers removed			
		Counter+=1
		self.ProximityToSurface(Volume,SurfaceCoords,ThresholdArray,Counter)


	def ExtractSurfaceCoordinates(self,Volume):
		#Get the Surface File
		Surface=self.GetSurface(Volume)

		#Get the Number of Points
		N=Surface.GetNumberOfPoints()
		
		#Loop over all of the points and extract the coordinates
		SurfaceCoordsX=np.zeros(N)
		SurfaceCoordsY=np.zeros(N)
		SurfaceCoordsZ=np.zeros(N)
		for i in range(N):
			surface_coords_=Surface.GetPoint(i)[:]
			SurfaceCoordsX[i]=surface_coords_[0]	
			SurfaceCoordsY[i]=surface_coords_[1]	
			SurfaceCoordsZ[i]=surface_coords_[2]	

		return [SurfaceCoordsX,SurfaceCoordsY,SurfaceCoordsZ]

        
	def Read_Vtu(self,filename):
		reader=vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(filename)
		reader.Update()
		return reader.GetOutput()   
        
	def WriteOutput(self,filename,data):
		Volume=vtk.vtkXMLUnstructuredGridWriter()
		Volume.SetFileName(filename)
		Volume.SetInputData(data)
		Volume.Update()
 
	def PRINT_PROGRESS(self,i,N,progress_old):
		progress_=(int((float(i)/N*100+0.5)))
		if progress_%10==0 and progress_%10!=progress_old: print ("    Progress: %d%%"%progress_)
		return progress_%10

	def GetSurface(self,Volume):
		Surface=vtk.vtkDataSetSurfaceFilter()
		Surface.SetInputData(Volume)
		Surface.Update()
		Surface=Surface.GetOutput()
		return Surface

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will scrap over the MBF vales from the interior and exterior edges.")

        #Input filename of the perfusion map
	parser.add_argument('-ifile', '--InputFileName', type=str, required=True, dest="InputFileName",help="Volumetric Mesh that contains the Myocardial Blood Flow Data and Territory Maps")

	#Option to strip off one endo- and one epi-cardial element
	parser.add_argument('-Layers', '--Layers', type=float, required=True, default=1, dest="Layers", help="The Number of Layers to Shed of the Myocardium")	

        #Array Name of the Data
	parser.add_argument('-arrayname', '--ArrayName', type=str, required=False,default="ImageScalars", dest="ArrayName",help="The name of the array containing the MBF values")

        #Output argumenets
	parser.add_argument('-ofile', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename to store the vtu file with boundary values scrapped off")

	
	args=parser.parse_args()

	ImageAnalysisMyocardiumRemoveBoundaryData(args).main()	
