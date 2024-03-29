#This script will read the territories file
#and write a tecplot file with histogram plots

import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
import argparse

class CoronaryPerfusionPlot():
	def __init__(self,Args):
                #Input argumenets
		self.Args=Args
		self.MBF_data={"RCA":[], "LAD":[], "Diag1":[], "Diag2":[], "Diag3":[], "LCx":[], "Ramus":[]}
	
		#The LV filename that contains MBF
		self.VentricleFileName=self.Args.InputFileName

                #Coronary Centerline Files
		CenterlineFileNamesPaths=sorted(glob("%s/wall_*cl.vtp"%self.Args.CenterlinesFolder))
		CenterlineFileNames=[filename.split("/")[-1] for filename in CenterlineFileNamesPaths]

                #Separate the Coronary Territories into Left and Right side 
		self.CenterlineFileNamesPaths=[]
		self.CenterlineFileNames=[]
		for i in range(len(CenterlineFileNames)):
			if CenterlineFileNames[i].find("RCA")>=0:
				self.CenterlineFileNames.append(CenterlineFileNames[i])
				self.CenterlineFileNamesPaths.append(CenterlineFileNamesPaths[i])
			if CenterlineFileNames[i].find("LCA")>=0:
				self.CenterlineFileNames.append(CenterlineFileNames[i])
				self.CenterlineFileNamesPaths.append(CenterlineFileNamesPaths[i])

	def main(self):
		#Read the LV mesh file that contains the territory mappings and MBF
		print ("--- Reading Left Ventricle Myocardial Blood Flow Data and Territory Maps: %s"%self.VentricleFileName)
		VentricleMesh=self.Read_Vtu(self.VentricleFileName)
		
		#Compute the vessel-specific MBF
		MBF_vessel=self.compute_mbf(VentricleMesh)
		
		#Write the tecplot file with the Vessel S 
		self.write_MBF_data(MBF_vessel)
	
	def compute_mbf(self,VentricleMesh):
                #Separate the vessels into various territories
		MBF_data=self.MBF_data
		
		#Get the number of points in the array
		N=VentricleMesh.GetNumberOfPoints()

                #Loop over all points and append to the correct array
		for i in range(N):
			value_=VentricleMesh.GetPointData().GetArray(self.Args.ArrayName).GetValue(i)
			territory_idx_=VentricleMesh.GetPointData().GetArray("TerritoryMaps").GetValue(i)
			filename_=self.CenterlineFileNames[territory_idx_]
			if filename_.find("RCA")>=0: MBF_data["RCA"].append(value_)
			elif filename_.find("lad")>=0: MBF_data["LAD"].append(value_)
			elif filename_.find("diag")>=0:MBF_data["Diag"].append(value_)
			elif filename_.find("septal")>=0: MBF_data["Septal"].append(value_)
			elif filename_.find("lcx")>=0: MBF_data["LCx"].append(value_)
			elif filename_.find("LCA")>=0: MBF_data["LAD"].append(value_)
			else:
				print ("A value was not found for any territory")
				exit(1)
			if filename_.find("LCA")>=0 and filename_.find("lcx")<0:
				MBF_data["LAD_Diag_Sep"].append(value_)
		return MBF_data

        
	def write_MBF_data(self,MBF_data):
                
		outfile=open(self.Args.OutputFileName,'w')
		outfile.write('TITLE="MBF(mL/min/100g)"\n')
		outfile.write('VARIABLES = "Region","MBF"\n')
                
		if self.Args.PostCABG==0: 
			count=1
			tag="pre"
		else: 
			count=1.25
			tag="post"

		for key,elem in MBF_data.items():
			outfile.write('Zone T= "%s_%s", I=3, F=POINT\n'%(key,tag))
			mean_=np.mean(elem)
			stdev_=np.std(elem)
			outfile.write("%.02f %.05f\n"%(count,mean_+stdev_))
			outfile.write("%.02f %.05f\n"%(count,mean_))
			outfile.write("%.02f %.05f\n"%(count,mean_-stdev_))
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
	parser = argparse.ArgumentParser(description="This script will generate tecplot files from perfusion territory vtu files.")

        #Input filename of the perfusion map
	parser.add_argument('-ifile', '--InputFileName', type=str, required=False, dest="InputFileName",help="Surface Mesh that contains the Myocardial Blood Flow Data and Territory Maps")

        #Input folder that contains the coronary centerlines
	parser.add_argument('-centerlinesfolder', '--CenterlinesFolder', type=str, required=True, dest="CenterlinesFolder",help="Folder that contains the centerline files")

        #Array Name of the Data
	parser.add_argument('-arrayname', '--ArrayName', type=str, required=False,default="ImageScalars", dest="ArrayName",help="The name of the array containing the MBF values")
	
	#Option to strip off one endo- and one epi-cardial element
	parser.add_argument('-stripedges', '--StripEdges', type=int, required=False, default=0, dest="StripEdges", help="The option will strip MBF values from exterior and interior of the myocardium")
	parser.add_argument('-edgespacing', '--EdgeSpacing', type=float, required=False, default=0, dest="EdgeSpacing", help="The distance from interior wall to remove MBF values")	

	#Pre- or Post-CABG
	parser.add_argument('-postCABG', '-PostCABG', type=int, required=True, default=0, dest="PostCABG", help="0 for pre-cabg (default) and 1 for post-cabg")

        #Output argumenets
	parser.add_argument('-ofile', '--OutputFileName', type=str, required=True, dest="OutputFileName",help="The output filename to store the perfusion plots and histograms")
	
	args=parser.parse_args()
	
	if args.StripEdges!=0 and args.EdgeSpacing==0:
		parser.error('EdgeSpacing not provided while StripeEdge is turned on')

	CoronaryPerfusionPlot(args).main()	
