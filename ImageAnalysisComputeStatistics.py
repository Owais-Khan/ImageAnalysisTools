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

class ImageAnalysisComputeStatistics():
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
		Npts=Volume.GetNumberOfPoints()
		DataArray=np.zeros(Npts)
		for i in range(Npts):
			DataArray[i]=Volume.GetPointData().GetArray(self.Args.ArrayName).GetValue(i)
			
		#Perform Statistics
		Mean=np.average(DataArray)
		Stdev=np.std(DataArray)
		Median50=np.percentile(DataArray,50)
		Median75=np.percentile(DataArray,75)
		Median80=np.percentile(DataArray,80)
		Median85=np.percentile(DataArray,85)
		Median90=np.percentile(DataArray,90)
		Median95=np.percentile(DataArray,95)

		print ("Writing Statistics: %s/MBF_Statistics.dat"%self.Args.OutputFolder)
		#Write the data to a file
		outfile=open("%s/MBF_Statistics.dat"%self.Args.OutputFolder,'w')
		outfile.write("Mean    MBF: %.05f\n"%Mean)
		outfile.write("Stdev   MBF: %.05f\n"%Stdev)
		outfile.write("Median  MBF: %.05f\n"%Median50)
		outfile.write("Perc75  MBF: %.05f\n"%Median75)
		outfile.write("Perc80  MBF: %.05f\n"%Median80)
		outfile.write("Perc85  MBF: %.05f\n"%Median85)
		outfile.write("Perc90  MBF: %.05f\n"%Median90)
		outfile.write("Perc95  MBF: %.05f\n"%Median95)
		outfile.write("------------------------------\n")
		outfile.write("MBF_norm_50: %.05f\n"%np.average(DataArray/Median50))
		outfile.write("MBF_norm_75: %.05f\n"%np.average(DataArray/Median75))
		outfile.write("MBF_norm_80: %.05f\n"%np.average(DataArray/Median80))
		outfile.write("MBF_norm_85: %.05f\n"%np.average(DataArray/Median85))
		outfile.write("MBF_norm_90: %.05f\n"%np.average(DataArray/Median90))
		outfile.write("MBF_norm_95: %.05f\n"%np.average(DataArray/Median95))
		outfile.close()

		print ("Writing Raw Data: %s/Data_Raw.dat"%self.Args.OutputFolder)
		outfile=open("%s/Data_Raw.txt"%self.Args.OutputFolder,'w')
		for i in range(Npts):
			outfile.write("%d %.05f\n"%(i,DataArray[i]))
		outfile.close()


		#Compute the number of Bins
		if self.Args.Bins is None:
			q25,q75=np.percentile(DataArray,[25,75])
			bin_width = int(2*(q75-q25)*len(DataArray)**(-1/3))
			self.Args.Bins=int(round((DataArray.max()-DataArray.min())/bin_width))
			print ("Freedman-Diaconis number of bins:",self.Args.Bins)
		
			
	
		plt.hist(DataArray,bins=self.Args.Bins,normed=True,label="Data",color='b',alpha=0.1)
		plt.xlabel("X")
		plt.ylabel("Probability")

		#Plot the Probibility Distribution Function
		if (self.Args.XMin is None) or (self.Args.XMax is None):
			self.Args.XMin,self.Args.XMax=plt.xlim()
			plt.xlim(XMin,XMax)
		else:
			plt.xlim(self.Args.XMin, self.Args.XMax)

		if (self.Args.YMin is None) or (self.Args.YMax is None):
			plt.ylim(0,0.02)
		else:
			plt.ylim(self.Args.YMin,self.Args.YMax)

		kde_xs=np.linspace(self.Args.XMin,self.Args.XMax,300)
		kde=st.gaussian_kde(DataArray)
		kde_pdf=kde.pdf(kde_xs)
		plt.plot(kde_xs,kde.pdf(kde_xs),label="PDF",color='k')
		

                #Plot the mean 
		plt.axvline(Mean,color='r',linestyle='-',label="Mean=%.01f\nStdev=%.01f"%(Mean,np.std(DataArray)))
                
                #Plot the mode
		Mode=kde_xs[np.argmax(kde_pdf)]
		plt.axvline(Mode,color='g',linestyle='-',label="Mode=%.01f"%Mode)


		plt.grid(linestyle = '--', linewidth = 0.5)

		plt.legend(loc="upper left")


		print ("Saving Figure to %s/Data_Figure.png"%self.Args.OutputFolder)	
		plt.savefig("%s/Data_Figure.png"%self.Args.OutputFolder,bbox_inches='tight',transparent=True)



if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute data statistics on MBF volume")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the myocardial volume")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName", default="scalars", help="The array name where the data is stored")
	
	parser.add_argument('-Bins', '--Bins', type=int, required=False, dest="Bins", help="The number of bins for the histogram.")
	parser.add_argument('-XMin', '--XMin', type=float, required=True, dest="XMin", help="The minimum for the plot.")
	parser.add_argument('-YMin', '--YMin', type=float, required=False, dest="YMin", help="The minimum for the plot.")
	parser.add_argument('-XMax', '--XMax', type=float, required=True, dest="XMax", help="The maximum for the plot.")
	parser.add_argument('-YMax', '--YMax', type=float, required=False, dest="YMax", help="The maximum for the plot.")

	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")

	

	args=parser.parse_args()
	ImageAnalysisComputeStatistics(args).Main()
