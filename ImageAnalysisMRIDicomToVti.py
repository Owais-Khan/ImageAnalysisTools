#!/Users/mokhan/Softwares/miniconda3/bin/python

#This file will create a velocity vtu file
#from phase and magnitude data from 4D MRI.
import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import numpy as np
import vtk
import argparse

class ImageAnalysisMRIDicomToVti():
	def __init__(self,Args):
                #Store the input arguments
		self.Args=Args

                #Load all of the phase file name
		filenames=glob("%s/*.dcm"%self.Args.InputFolderName)
		filenames=glob("%s/*.IMA"%self.Args.InputFolderName)

		#Sort all of the filenames into order from phase1-3 and magnitude
		#Find the first integer value in the filename
		filename_=filenames[0].split(".")
		for i in range(len(filename_)): 
			if filename_[i].isnumeric():
				sorttag=i 
				break
		filenames=sorted(filenames, key=lambda filename: int(filename.split(".")[sorttag]))
		print ("There are %d files in the folder"%len(filenames))
		
		#Now separate each series into seperate bins
		N_file=int(len(filenames)/4)
		print ("Sperating the files into %d chucks for phase1,2,3 and magnitude"%N_file)
		filenames_dic={"Phase1":[],"Phase2":[],"Phase3":[],"Magnitude":[]}
		filenames_dic["Phase1"]   =filenames[0:N_file]
		filenames_dic["Phase2"]   =filenames[N_file:2*N_file]
		filenames_dic["Phase3"]   =filenames[2*N_file:3*N_file]
		filenames_dic["Magnitude"]=filenames[3*N_file:]

		print ("Sorting each phase and magnitude from t=0 to t=N and then slice #")	
		#Now sort each bin from t=0 to t=N.
		filenames_dic["Phase1"]=sorted(filenames_dic["Phase1"],key=lambda filename: int(filename.split(".")[-2]))
		filenames_dic["Phase2"]=sorted(filenames_dic["Phase2"],key=lambda filename: int(filename.split(".")[-2]))
		filenames_dic["Phase3"]=sorted(filenames_dic["Phase3"],key=lambda filename: int(filename.split(".")[-2]))
		filenames_dic["Magnitude"]=sorted(filenames_dic["Magnitude"],key=lambda filename: int(filename.split(".")[-2]))
		
		self.FileNames=filenames_dic		
	
	def main(self):
		#Loop over all of the phases and magnitudes
		os.system("mkdir %s/temp_folder/"%self.Args.OutputFolderName)
		for key in self.FileNames:
			print ("--- Looping over:%s"%key)
			#Chunck up into N timesteps for each slice
			filenames_=list(self.Chunks(self.FileNames[key],self.Args.NumberOfTimesteps))
			
			#Now transpose so we have each timestep domain seperate	
			filenames_=list(map(list, zip(*filenames_)))
			
			#Create a shortcut link in temp directory to allow DicomReader a folder
			for i in range(self.Args.NumberOfTimesteps):
				#empty the temp_folder
				os.system("rm %s/temp_folder/*"%self.Args.OutputFolderName)
				for filename_ in filenames_[i]:
					filename_strip_=filename_.split("/")[-1].replace(".IMA",".dcm")
					os.system("cp %s %s/temp_folder/%s"%(filename_,self.Args.OutputFolderName,filename_strip_))
				#Write a vti file for each timestep 				
				self.ConvertDicomDirectoryToVti(self.Args.OutputFolderName+"/temp_folder","%s/%s_timestep%.04d.vti"%(self.Args.OutputFolderName,key,i))

	def ConvertDicomDirectoryToVti(self,inputfoldername,outfilename):
		print ("---------- Writing %s"%outfilename)
		#Read the dicom files in the directory
		reader=vtk.vtkDICOMImageReader()
		reader.SetDirectoryName(inputfoldername)
		reader.Update()
		#Write the file in vti format
		writer=vtk.vtkXMLImageDataWriter()
		writer.SetFileName(outfilename)
		writer.SetInputData(reader.GetOutput())
		writer.Write()


	def Chunks(self,lst, n):
		"""Yield successive n-sized chunks from lst."""
		for i in range(0, len(lst), n):
			yield lst[i:i + n]
	
					
				

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take a dicom folder with phase 1,2,3 and write a velocity vector file")

        #Input filename of the Phase maps
	parser.add_argument('-ifolder', '--InputFolderName', type=str, required=True, dest="InputFolderName",help="The name of the folder contain three phases and magnitude")
        
        #Number of Timesteps per cycle
	parser.add_argument('-timesteps', '--NumberOfTimesteps', type=int, required=True, dest="NumberOfTimesteps",help="The number of time points stored per cardiac cycle")
        
	#The length/period of cardiac cycle
	parser.add_argument('-period', '--Period', type=float, required=False, default=1, dest="Period",help="The length/period of the cardiac cycle")

        #Diastolic cycle
	parser.add_argument('-venc', '--VelocityEncoding', type=float, required=True, dest="VelocityEncoding",help="The value of velocity encoding used for the MRI scan")

        #Output folder
	parser.add_argument('-ofolder', '--OutputFolderName', type=str, required=True, dest="OutputFolderName",help="The name of the folder to store the velocity data")

	args=parser.parse_args()
	ImageAnalysisMRIDicomToVti(args).main()

