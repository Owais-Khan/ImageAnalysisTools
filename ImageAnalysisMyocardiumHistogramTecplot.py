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
			self.Args.OutputFolder=self.Args.InputFileName.replace(self.Args.InputFileName.split("/")[-1],"/")

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
		DataArray=vtk_to_numpy(Volume.GetPointData().GetArray(self.Args.ArrayName)).astype(float)
		
		#Compute Mean and standard deviation
		print ("--- Computing Mean and Standard Deviation")
		mu, std = st.norm.fit (DataArray)
		print ("--- Computing 75th Percentile")
		MBF75Perc=np.percentile(DataArray,75)
		print ("--- Computing the Median")
		MBF50Perc=np.percentile(DataArray,50)
		print ("--- Performing Shapiro-Wilks test for Normality")
		ShapiroPvalue=st.shapiro(DataArray)[1]
		print ("--- Performing Kolmogorov Smirnov test for Normality")
		KolmogorovPvalue=st.kstest(DataArray,'norm')[1]
		print ("--- Computing Kurtosis and Skewness")
		Kurtosis=st.kurtosis(DataArray)
		Skewness=st.skew(DataArray)
		
		#Compute the region of Left ventricle below 100ml/g/min (MBF) and 0.8 (IndexMBF)
		counter=0.0
		MBF_counter=0
		IndexMBF_counter=0
		for i in range(len(DataArray)):
			if DataArray[i]==0.0: continue
			else:
				if DataArray[i]<100: MBF_counter+=1
				if (DataArray[i]/MBF75Perc)<0.8: IndexMBF_counter+=1
				counter+=1.0

		LV_Below_100_MBF=(float(MBF_counter)/counter)*100
		LV_Below_08_IndexMBF=(float(IndexMBF_counter)/counter)*100
			
		
		print ("--- Writing the Statistics to An OutputFile")
		outfile=open("%s/MBF_statistics.txt"%self.Args.OutputFolder,'w')
		outfile.write("Mean MBF:                    %.05f\n"%mu)	
		outfile.write("Std MBF:                     %.05f\n"%std)	
		outfile.write("Mean Index MBF:              %.05f\n"%(mu/MBF75Perc))	
		outfile.write("Std Index MBF:               %.05f\n"%(std/MBF75Perc))	
		outfile.write("50Percentile:                %.05f\n"%MBF50Perc)	
		outfile.write("75Percentile:                %.05f\n"%MBF75Perc)	
		outfile.write("ShapriWilks p-value:         %.05f\n"%ShapiroPvalue)
		outfile.write("Kolmogorov Smirnov p-value:  %.05f\n"%KolmogorovPvalue)
		outfile.write("Kurtosis:                    %.05f\n"%Kurtosis)
		outfile.write("Skewness:                    %.05f\n"%Skewness)
		outfile.write("LV Below 100ml/g/min MBF:    %.05f%%\n"%LV_Below_100_MBF)
		outfile.write("LV Below 0.8 Index MBF:      %.05f%%\n"%LV_Below_08_IndexMBF)
		outfile.close()
	
		#Compute Histogram
		x_axis=np.linspace(DataArray.min(),DataArray.max(),self.Args.Bins)
		kde=st.gaussian_kde(DataArray)
		kde_pdf=kde.pdf(x_axis)

		#Now write the Tecplot Verson of the data
		print ("--- Writing Tecplot Histogram in %s/MBF_HistogramTecplot.dat"%self.Args.OutputFolder)
		outfile=open("%s/MBF_Histogram_Tecplot.dat"%self.Args.OutputFolder,'w')
		outfile.write('TITLE="MBF(mL/min/100g)"\n')
		outfile.write('VARIABLES = "Bins","MBF"\n')
		
		#First Write Normal Distribution of the Data
		outfile.write('Zone T= "NormalDistribution", I=%d, F=POINT\n'%self.Args.Bins)
		pdf_normal = st.norm.pdf(x_axis, mu, std)	
		for i in range(self.Args.Bins):
			outfile.write("%.05f %.05f\n"%(x_axis[i],pdf_normal[i]))

		#Now write the raw data
		outfile.write('Zone T= "DataDistribution", I=%d, F=POINT\n'%self.Args.Bins)
		for i in range(self.Args.Bins):
			outfile.write("%.05f %.05f\n"%(x_axis[i],kde_pdf[i]))
			
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
