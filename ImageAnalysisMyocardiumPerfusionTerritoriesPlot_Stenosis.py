#This script will read the territories file
#and write a tecplot file with histogram plots

import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
import argparse

class ImageAnalysisMyocardiumPerfusionTerritoriesPlot():
	def __init__(self,Args):
                #Input argumenets
		self.Args=Args

                #If output filename is not defined, create your own
		if self.Args.OutputFileName is None:
			InputFileNameStripped=self.Args.InputFileName.split("/")[-1]
			OutputFolder =self.Args.InputFileName.replace(InputFileNameStripped,"")
			OutputFileName=OutputFolder+"MyocardiumPerfusionTerritoriesPlot_Stenosis.dat"
			self.Args.OutputFileName=OutputFileName

		#The LV filename that contains MBF
		self.VentricleFileName=self.Args.InputFileName

                #Coronary Centerline Files
		self.CenterlineFileNamesPaths=sorted(glob("%s/*_Stenosis_*.vtp"%self.Args.CenterlinesFolder))
		self.CenterlineFileNames=[filename.split("/")[-1] for filename in self.CenterlineFileNamesPaths]
			

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
		MBF_data={}
		
		#Get the number of points in the array
		N=VentricleMesh.GetNumberOfPoints()

                #Loop over all points and append to the correct array
		for i in range(N):
			value_=VentricleMesh.GetPointData().GetArray(self.Args.ArrayName).GetValue(i)
			territory_idx_=VentricleMesh.GetPointData().GetArray("TerritoryMaps").GetValue(i)
			filename_=self.CenterlineFileNames[territory_idx_]
			CoronaryName_=filename_.split("_")[0]
			if CoronaryName_ in MBF_data:
				MBF_data[CoronaryName_].append(value_)
			else: 
				MBF_data[CoronaryName_]=[]
				MBF_data[CoronaryName_].append(value_)

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
		reader=vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(filename)
		reader.Update()
		return reader.GetOutput()

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will generate tecplot files from perfusion territory vtu files.")

        #Input filename of the perfusion map
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="Volumetric Mesh that contains the Myocardial Blood Flow Data and Territory Maps")

        #Input folder that contains the coronary centerlines
	parser.add_argument('-CenterlinesFolder', '--CenterlinesFolder', type=str, required=True, dest="CenterlinesFolder",help="Folder that contains the centerline files")

        #Array Name of the Data
	parser.add_argument('-arrayname', '--ArrayName', type=str, required=False,default="scalars", dest="ArrayName",help="The name of the array containing the MBF values")

	#Pre- or Post-CABG
	parser.add_argument('-postCABG', '-PostCABG', type=int, required=False, default=0, dest="PostCABG", help="0 for pre-cabg (default) and 1 for post-cabg")

        #Output argumenets
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename to store the perfusion plots and histograms")

	args=parser.parse_args()
	ImageAnalysisMyocardiumPerfusionTerritoriesPlot(args).main()	
