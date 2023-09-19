import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *
import matplotlib.pyplot as plt
import csv
import scipy.stats as st

class ImageAnalysisHistogramTecplot():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			os.system("mkdir Results")
			self.Args.OutputFolder="./Results/"

	def Main(self):
                #Read the vtu file
		print ("--- Reading %s"%self.Args.InputFileName)
		if   self.Args.InputFileName[-4:]==".vtk":
			Volume=ReadVTKFile(self.Args.InputFileName)   #Read a VTK volume stack
			Volume=ThresholdByUpper(Volume,self.Args.ArrayName,1) #Convert to an unstructured grid 
                
		elif self.Args.InputFileName[-4:]==".vtu":
			Volume=ReadVTUFile(self.Args.InputFileName) #Read a VTU unstructured grid
		else:
			print ("The extension %s is not valid for volume"%self.Args.InputFileName[-4:])
			print ("Exiting...")
			exit(1)

		
		#Compute All of the Scalar quantities
		print ("--- Reading Scalar Values")
		DataArray=vtk_to_numpy(Volume.GetPointData().GetArray(self.Args.ArrayName)
		
		#Compute Mean and standard deviation
		mu, std = st.norm.fit (DataArray)
	
		#Compute Histogram
		x_axis=np.linspace(DataArray.min(),DataArray.max(),self.Args.Bins)
		kde=st.gaussian_kde(DataArray)
		kde_pdf=kde.pdf(x_axis)

		#Now write the Tecplot Verson of the data
		outfile=open("%s/MBF_Statistics_Tecplot.dat"%self.Args.OutputFolder)
		outfile.write('TITLE="MBF(mL/min/100g)"\n')
		outfile.write('VARIABLES = "Bins","MBF"\n')
		
		#First Write Normal Distribution of the Data
		outfile.write('Zone T= "NormalDistribution", I=%d, F=POINT\n'%self.Args.Bins)
		pdf_normal = st.norm.pdf(x_axis, mu, std)	
		for i in range(self.Args.Bins):
			outfile.write("%.05f %.05f\n"%(x_axis[i],pdf_normal[i])

		#Now write the raw data
		outfile.write('Zone T= "DataDistribution", I=%d, F=POINT\n'%self.Args.Bins)
		for i in range(self.Args.Bins):
			outfile.write("%.05f %.05f\n"%(x_axis[i],kde_pdf[i])
			
		#Now write the Mean and Standard Deviation Line
		outfile.write('Zone T= "Mean", I=2, F=POINT\n')
		outfile.write('%.05f %.05f\n'%(mu,0))
		outfile.write('%.05f %.05f\n'%(mu,1000))
		
		outfile.write('Zone T= "Std_negative", I=2, F=POINT\n')
		outfile.write('%.05f %.05f\n'%(mu-std,0))
		outfile.write('%.05f %.05f\n'%(mu-std,1000))
		
		outfile.write('Zone T= "Std_negative", I=2, F=POINT\n')
		outfile.write('%.05f %.05f\n'%(mu+std,0))
		outfile.write('%.05f %.05f\n'%(mu+std,1000))

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute data statistics on MBF volume")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the myocardial volume")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName", default="scalars", help="The array name where the data is stored")
	
	parser.add_argument('-Bins', '--Bins', type=int, required=False, default=300, dest="Bins", help="The number of bins for the histogram.")

	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")

	

	args=parser.parse_args()
	ImageAnalysisHistogramTecplot(args).Main()
