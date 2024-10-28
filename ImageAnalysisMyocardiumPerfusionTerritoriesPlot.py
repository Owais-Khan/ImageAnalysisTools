#This script will read the territories file
#and write a tecplot file with histogram plots

import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
from scipy.stats import iqr as IQR
import argparse

class ImageAnalysisMyocardiumPerfusionTerritoriesPlot():
	def __init__(self,Args):
                #Input argumenets
		self.Args=Args

                #If output filename is not defined, create your own
		if self.Args.OutputFileName is None:
			InputFileNameStripped=self.Args.InputFileName.split("/")[-1]
			OutputFolder =self.Args.InputFileName.replace(InputFileNameStripped,"")
			OutputFileName=OutputFolder+"MyocardiumPerfusionTerritoriesPlot.dat"
			OutputFileName2=OutputFolder+"MyocardiumPerfusionTerritories_Labels.dat"
			self.Args.OutputFileName=OutputFileName

		#Name of the Coronary Artery Tags
		#self.MBF_data={"RCA":[], "LAD":[], "Diag":[], "Septal":[], "LCx":[], "LAD_Diag_Sep":[]}
		self.MBF_data={"LAD":[], "LCx":[], "Intermedius":[], "Diag1":[], "Diag2":[], "PDA":[], "PLA":[]}
		self.MBF_Labels={"LAD":[], "LCx":[], "Intermedius":[], "Diag1":[], "Diag2":[], "PDA":[], "PLA":[]}

		
		#Read the Territory Labels
		infile=open(self.Args.InputFileName.replace(".vtu","_Labels.dat"),'r')
		infile.readline()
		for LINE in infile:
			line=LINE.split()
			for key_ in self.MBF_Labels.keys():
				if line[1].find(key_)>=0: self.MBF_Labels[key_].append(int(line[0]))

		#The LV filename that contains MBF
		self.VentricleFileName=self.Args.InputFileName

                #Find the array name with work "scalars"
		if self.Args.ArrayName is None:
			for i in range(VentricleMesh.GetPointData().GetNumberOfArrays()):
				ArrayName_=str(VentricleMesh.GetPointData().GetArrayName(i))
				if ArrayName_.find("calars")>0: self.Args.ArrayName=ArrayName_
		print ("--- The array that contains MBF values is: %s"%self.Args.ArrayName)


	def main(self):
		#Read the LV mesh file that contains the territory mappings and MBF
		print ("--- Reading Left Ventricle Myocardial Blood Flow Data and Territory Maps: %s"%self.VentricleFileName)
		VentricleMesh=self.Read_Vtu(self.VentricleFileName)
		
		#Compute the vessel-specific MBF
		MBF_vessel=self.compute_mbf(VentricleMesh)

		#Compute 75th Percentile for Normalization
		MBF_Data_LV=vtk_to_numpy(VentricleMesh.GetPointData().GetArray(self.Args.ArrayName))
		MeanLV=np.mean(MBF_Data_LV)
		StdevLV=np.std(MBF_Data_LV)
		MedianLV=np.median(MBF_Data_LV)
		IqrLV=IQR(MBF_Data_LV)
		Per25LV=np.percentile(MBF_Data_LV,25)
		Per75LV=np.percentile(MBF_Data_LV,75)
		MaxLV=np.max(MBF_Data_LV)
		MinLV=np.min(MBF_Data_LV)

		
		#Write the tecplot file with the Vessel S 
		self.write_MBF_data(MBF_vessel)

		#Writing the Data to ASCII format
		print ("--- Writing the Territory Statistics: %s"%self.Args.OutputFileName.replace(".dat","_statistics.dat"))
		outfile=open(self.Args.OutputFileName.replace(".dat","_statistics.dat"),'w')
                                
		outfile.write("Entire Myocardium\n")
		outfile.write("---Mean:          %.06f\n"%MeanLV)
		outfile.write("---Stdev:         %.06f\n"%StdevLV)
		outfile.write("---Median:        %.06f\n"%MedianLV)
		outfile.write("---IQR:           %.06f\n"%IqrLV)
		outfile.write("---25th Percentile: %.06f\n"%Per25LV)
		outfile.write("---75th Percentile: %.06f\n"%Per75LV)
		outfile.write("---Max:           %.06f\n"%MaxLV)
		outfile.write("---Min:           %.06f\n"%MinLV)
		outfile.write("\n")


		for key_ in self.MBF_data.keys():
			if len(self.MBF_Labels[key_])>0:
				DataArray_=np.array(self.MBF_data[key_])
				outfile.write("%s\n"%key_)
				outfile.write("---Mean:            %.06f\n"%np.mean(self.MBF_data[key_]))
				outfile.write("---Stdev:           %.06f\n"%np.std(self.MBF_data[key_]))
				outfile.write("---MeanNormalized:  %.06f\n"%(np.mean(self.MBF_data[key_])/Per75LV))
				outfile.write("---StdevNormalized: %.06f\n"%(np.std(self.MBF_data[key_])/Per75LV))
				outfile.write("---Median:          %.06f\n"%np.median(self.MBF_data[key_]))
				outfile.write("---IQR:             %.06f\n"%IQR(self.MBF_data[key_]))
				outfile.write("---25th Percentile: %.06f\n"%np.percentile(self.MBF_data[key_],25))
				outfile.write("---75th Percentile: %.06f\n"%np.percentile(self.MBF_data[key_],75))
#				outfile.write("---Mode:    %.06f\n"%sp.stats.mode(self.MBF_data[key_]))
				outfile.write("---Max:             %.06f\n"%np.max(self.MBF_data[key_]))
				outfile.write("---Min:             %.06f\n"%np.min(self.MBF_data[key_]))

				#Compute the percent volume of the territory
				VolumePercent=len(self.MBF_data[key_])/VentricleMesh.GetNumberOfPoints()*100
				outfile.write("---%%Volume: %.06f\n"%VolumePercent)
				outfile.write("\n")
	
	def compute_mbf(self,VentricleMesh):
                #Separate the vessels into various territories
		MBF_data=self.MBF_data
		
		#Get the number of points in the array
		N=VentricleMesh.GetNumberOfPoints()

                #Loop over all points and append to the correct array
		for i in range(N):
			value_=VentricleMesh.GetPointData().GetArray(self.Args.ArrayName).GetValue(i)
			if value_==0: continue
			territory_idx_=VentricleMesh.GetPointData().GetArray("TerritoryMaps").GetValue(i)
			for key_ in self.MBF_Labels.keys():
				if int(territory_idx_) in self.MBF_Labels[key_]: self.MBF_data[key_].append(value_)	
			
		return MBF_data


        
	def write_MBF_data(self,MBF_data):
		print ("--- Writing the Data for Tecplot: %s"%self.Args.OutputFileName)                
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
			if len(self.MBF_Labels[key])>0:
				mean_=np.mean(elem)
				stdev_=np.std(elem)
			else:
				mean_=0
				stdev_=0
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

        #Array Name of the Data
	parser.add_argument('-Arrayname', '--ArrayName', type=str, required=False,default="ImageScalars", dest="ArrayName",help="The name of the array containing the MBF values")

	#Pre- or Post-CABG
	parser.add_argument('-PostCABG', '-PostCABG', type=int, required=False, default=0, dest="PostCABG", help="0 for pre-cabg (default) and 1 for post-cabg")

        #Output argumenets
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename to store the perfusion plots and histograms")

	args=parser.parse_args()
	ImageAnalysisMyocardiumPerfusionTerritoriesPlot(args).main()	
