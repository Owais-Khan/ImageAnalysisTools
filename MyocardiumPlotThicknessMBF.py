import numpy as np
import os
import sys
import argparse
import matplotlib.cm as cm
import matplotlib.pyplot as plt

class MyocardiumPlotThicknessMBF():
	def __init__(self,args):
		self.Args=args
		if self.Args.OutputFile is None:
			self.Args.OutputFile=self.Args.InputFile.replace(".dat","_tecplot.dat")
		self.LengthRes=50
		self.ThicknessRes=10 
		self.Thickness=2.5 #Normalize the Thickness by 2cm
		self.Length=1. #The normalize length of ApexBase
	def Main(self):
		#Create a 2D Array of Apex-Base location and MBF Along the Thickness
		DataArray=np.zeros(shape=(self.LengthRes,self.ThicknessRes))	
		Counter =np.zeros(shape=(self.LengthRes,self.ThicknessRes))	
		
		#Read the data and assign to an array
		infile=open(self.Args.InputFile,'r')
		infile.readline()
		for LINE in infile:
			line_=LINE.split()
			#Get the Location on Apex-Base axis
			Yloc_=int(float(line_[0])*self.LengthRes+0.5)
			#Get the Number of Points along the wall thickness
			Npts_=int(line_[1])
			#Fill out the Number of points
			for i in range(0,Npts_,2):
				Xloc_=float(line_[i+2])/self.Thickness #Thickness normalized
				Xloc_=int((Xloc_*self.ThicknessRes+0.5))
				DataArray[Yloc_-1,Xloc_]+=float(line_[i+3])
				Counter[Yloc_-1,Xloc_]+=float(1)
		infile.close()
		
		#Fix the counter (i.e., add ones so averaging doesn't lead to 0/0) 	
		for i in range(0,self.LengthRes):
			for j in range(0,int(self.ThicknessRes)):
				if Counter[i,j]==0: Counter[i,j]=1		

		#Now Divide the DataArray by Counter Array
		DataArray=DataArray/Counter	

		self.WriteTecplotContour(DataArray)

	def WriteTecplotContour(self,DataArray):	
		#Write the Tecplot File
		outfile=open(self.Args.OutputFile,'w')
		outfile.write('VARIABLES= "ApexBaseLine", "Thickness", "MBF"\n')
		outfile.write('ZONE I=%d, J=%d, DATAPACKING=POINT\n'%(self.ThicknessRes,self.LengthRes))
		for i in range(0,self.LengthRes):
			I_loc_=(float(i)/self.LengthRes)*self.Length
			for j in range(0,self.ThicknessRes):
				J_loc_=(float(j))/self.ThicknessRes*self.Thickness
				outfile.write("%d %d %.05f\n"%(i,j,DataArray[i,j]))
		outfile.close()
		
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will plot the thickness values output my MyocardiumProjectVolumeToSurface. It will generate a tecplot file that can be load")

	parser.add_argument('-InputFile', '--InputFile', type=str, required=True, dest="InputFile",help="The text file that contains the thickness values")

        #Output Filename 
	parser.add_argument('-OutputFile', '--OutputFile', type=str, required=False, dest="OutputFile",help="The output file that contains the MBF thickness values")

        
	args=parser.parse_args()
	MyocardiumPlotThicknessMBF(args).Main()


