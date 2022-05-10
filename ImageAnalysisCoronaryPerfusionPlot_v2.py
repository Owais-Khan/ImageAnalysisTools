#This script will read the territories file
#and write a tecplot file with histogram plots

import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
import argparse
from utilities import *
import copy

class CoronaryPerfusionPlot():
	def __init__(self,Args):
                #Input argumenets
		self.Args=Args

		if self.Args.OutputFileName is None:
			self.Args.OutputFileName="Results/MBF_ProjectedToSurface_Territories_Plots.dat"

		self.Args.Array1="MBF_Normalized75Q"
		self.Args.Array2="MBF_Normalized50Q"
		self.Args.Array3="MBF_WallAveraged"
		self.Args.Array4="MBF_NormalizedAVG"
	
		self.MBF_data={"RCA":[], "LAD":[], "LCx":[], "Ramus":[], "Diag1":[] }
	
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
				print ("------ Loading Centerline: %s"%CenterlineFileNamesPaths[i])

	def Main(self):
		#Read the LV mesh file that contains the territory mappings and MBF
		print ("--- Reading Left Ventricle Myocardial Blood Flow Data and Territory Maps: %s"%self.VentricleFileName)
		VentricleMesh=ReadVTPFile(self.VentricleFileName)
		
		#Compute the vessel-specific MBF
		MBF_Array1,MBF_Array2,MBF_Array3,MBF_Array4=self.compute_mbf(VentricleMesh)
		
		#Write the tecplot file with the Vessel S 
		self.write_MBF_data(MBF_Array1,MBF_Array2,MBF_Array3,MBF_Array4)
	
	def compute_mbf(self,VentricleMesh):
                #Separate the vessels into various territories
		MBF_Array1=copy.deepcopy(self.MBF_data)
		MBF_Array2=copy.deepcopy(self.MBF_data)
		MBF_Array3=copy.deepcopy(self.MBF_data)
		MBF_Array4=copy.deepcopy(self.MBF_data)
			
		#Get the number of points in the array
		N=VentricleMesh.GetNumberOfPoints()

                #Loop over all points and append to the correct array
		for i in range(N):
			ValueArray1_=VentricleMesh.GetPointData().GetArray(self.Args.Array1).GetValue(i)
			ValueArray2_=VentricleMesh.GetPointData().GetArray(self.Args.Array2).GetValue(i)
			ValueArray3_=VentricleMesh.GetPointData().GetArray(self.Args.Array3).GetValue(i)
			ValueArray4_=VentricleMesh.GetPointData().GetArray(self.Args.Array3).GetValue(i)
	
			territory_idx_=VentricleMesh.GetPointData().GetArray("TerritoryMaps").GetValue(i)
			filename_=self.CenterlineFileNames[territory_idx_]
			
			for Key_ in self.MBF_data:
				if filename_.find(Key_)>=0:
					MBF_Array1[Key_].append(ValueArray1_)
					MBF_Array2[Key_].append(ValueArray2_)
					MBF_Array3[Key_].append(ValueArray3_)
					MBF_Array4[Key_].append(ValueArray4_)
			#else:
				#print ("A value was not found for %s"%Key_)
				#print ("Removing %s from the keys")
				#exit(1)
		return MBF_Array1,MBF_Array2,MBF_Array3,MBF_Array4

        
	def write_MBF_data(self,MBF_Array1,MBF_Array2,MBF_Array3,MBF_Array4):
                
		outfile=open(self.Args.OutputFileName,'w')
		outfile.write('TITLE="MBF(mL/min/100g)"\n')
		outfile.write('VARIABLES = "Region","%s", "%s", "%s", "%s"\n'%(self.Args.Array1,self.Args.Array2,self.Args.Array3,self.Args.Array4))
                
		if self.Args.PostCABG==0: 
			count=1
			tag="pre"
		else: 
			count=1.25
			tag="post"

		for key in self.MBF_data:
			outfile.write('Zone T= "%s_%s", I=3, F=POINT\n'%(key,tag))
			mean_Array1=np.mean(MBF_Array1[key])
			mean_Array2=np.mean(MBF_Array2[key])
			mean_Array3=np.mean(MBF_Array3[key])
			mean_Array4=np.mean(MBF_Array4[key])

			stdev_Array1=np.std(MBF_Array1[key])
			stdev_Array2=np.std(MBF_Array2[key])
			stdev_Array3=np.std(MBF_Array3[key])
			stdev_Array4=np.std(MBF_Array4[key])

			outfile.write("%.02f %.05f %.05f %.05f %.05f\n"%(count,mean_Array1+stdev_Array1,mean_Array2+stdev_Array2,mean_Array3+stdev_Array3,mean_Array4+stdev_Array4))
			outfile.write("%.02f %.05f %.05f %.05f %.05f\n"%(count,mean_Array1,mean_Array2,mean_Array3,mean_Array4))
			outfile.write("%.02f %.05f %.05f %.05f %.05f\n"%(count,mean_Array1-stdev_Array1,mean_Array2-stdev_Array2,mean_Array3-stdev_Array3,mean_Array4-stdev_Array4))
			count+=1
		outfile.close()

	def Read_Vtu(self,filename):
		Volume=vtk.vtkXMLUnstructuredGridReader()
		Volume.SetFileName(filename)
		Volume.Update()
		Volume=Volume.GetOutput()
		N_p_vol=Volume.GetNumberOfPoints()

		#This options strips a single node from the surface
		#to avoid artifacts close to the surface of endo- or
		#epicardium
		if self.Args.StripEdges!=0:
			#Get the surface from the volume
			Surface=vtk.vtkDataSetSurfaceFilter()
			Surface.SetInputData(Volume)
			Surface.Update()
			Surface=Surface.GetOutput()
			
			#Get all of the surface nodes:
			N_p_surf=Surface.GetNumberOfPoints()
			Surface_Coords_x=np.zeros(N_p_surf)
			Surface_Coords_y=np.zeros(N_p_surf)
			Surface_Coords_z=np.zeros(N_p_surf)
			for i in range(N_p_surf):
				surface_coord_=Surface.GetPoint(i)[:]
				Surface_Coords_x[i]=surface_coord_[0]
				Surface_Coords_y[i]=surface_coord_[1]
				Surface_Coords_z[i]=surface_coord_[2]
			
			Surface_Coords_list_x=np.ndarray.tolist(Surface_Coords_x)
			Surface_Coords_list_y=np.ndarray.tolist(Surface_Coords_y)
			Surface_Coords_list_z=np.ndarray.tolist(Surface_Coords_z)

			#Tag all of the surface coordinates on the volume
			#Make the values on the Surface equal to 0
			Volume_Surface_tags=[]
			for i in range(0,N_p_vol):
				volume_coord_=Volume.GetPoint(i)[:]
				if (volume_coord_[0] in Surface_Coords_list_x) and (volume_coord_[1] in Surface_Coords_list_y) and (volume_coord_[2] in Surface_Coords_list_z):
					Volume_Surface_tags.append(True)
					Surface_Coords_list_x.remove(volume_coord_[0])
					Surface_Coords_list_y.remove(volume_coord_[1])
					Surface_Coords_list_z.remove(volume_coord_[2])
					Volume_Surface_tags.append(True)
				else:
					Volume_Surface_tags.append(False)
			

			exit(1)
		return reader

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will generate tecplot files from perfusion territory surface files.")

        #Input filename of the perfusion map
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=False, dest="InputFileName", default="Results/MBF_ProjectedToSurface_Territories.vtp",help="Surface Mesh that contains the Myocardial Blood Flow Data and Territory Maps")

        #Input folder that contains the coronary centerlines
	parser.add_argument('-CenterlinesFolder', '--CenterlinesFolder', type=str, required=False, default="Results/Centerlines", dest="CenterlinesFolder",help="Folder that contains the centerline files")

	#Option to strip off one endo- and one epi-cardial element
	parser.add_argument('-stripedges', '--StripEdges', type=int, required=False, default=0, dest="StripEdges", help="The option will strip MBF values from exterior and interior of the myocardium")
	parser.add_argument('-edgespacing', '--EdgeSpacing', type=float, required=False, default=0, dest="EdgeSpacing", help="The distance from interior wall to remove MBF values")	

	#Pre- or Post-CABG
	parser.add_argument('-postCABG', '-PostCABG', type=int, required=True, default=0, dest="PostCABG", help="0 for pre-cabg (default) and 1 for post-cabg")

        #Output argumenets
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename to store the perfusion plots and histograms")
	
	args=parser.parse_args()
	
	if args.StripEdges!=0 and args.EdgeSpacing==0:
		parser.error('EdgeSpacing not provided while StripeEdge is turned on')

	CoronaryPerfusionPlot(args).Main()	
