import vtk
import numpy
import argparse
from glob import glob
import os

class ImageAnalysis4DMRIAngiography():
	def __init__(self,Arg):
		self.Args=Arg

	def Main(self):
		#Separate all images into individual folders
		self.SeparateCycles(self.Args.Phase1Folder,"Phase1")		

	def SeparateCycles(self,FolderName,Tag):
		#Get all of the files in the folder
		FileNames_=sorted(glob(FolderName+"/*.IMA"))+sorted(glob(FolderName+"/*.dcm"))
		
		#Separte the files in N folders
		os.system("mkdir %s/%s"%(self.Args.OutputFolder,Tag))
		for i in range(self.Args.Cycles): os.system("mkdir %s/%s/%s_%d"%(self.Args.OutputFolder,Tag,Tag,i))	
		#Loop over all of the files and stack them into difference folders
		counter=0
		for i in range(len(FileNames_)):
			os.system("cp %s %s/%s/%s_%d/"%(FileNames_[i],self.Args.OutputFolder,Tag,Tag,counter))
			if i%(len(FileNames_)/self.Args.Cycles): counter+=1
			




if __name__=="__main__":
	#Description
	parser = argparse.ArgumentParser(description="This script will average the 4D Flow MRI data to convert into angiography data to segmentation easier.")

	#Provide a path to the Magnitude Images
	parser.add_argument('-MagnitudeFolder', '--MagnitudeFolder', type=str, required=True, dest="MagnitudeFolder",help="The foldername that contains the magnitude images")
	
	#Provide a path for phase images 1 2 and 3
	parser.add_argument('-Phase1Folder', '--Phase1Folder', type=str, required=True, dest="Phase1Folder",help="The foldername that contains the phase1 images")
	parser.add_argument('-Phase2Folder', '--Phase2Folder', type=str, required=True, dest="Phase2Folder",help="The foldername that contains the phase2 images")
	parser.add_argument('-Phase3Folder', '--Phase3Folder', type=str, required=True, dest="Phase3Folder",help="The foldername that contains the phase3 images")
	
	#Provide a velocity encoding 
	parser.add_argument('-Venc', '--Venc', type=int, required=True, dest="Venc",help="The velocity encoding used for this scan.")

	#Number of scan cycles 
	parser.add_argument('-Cycles', '--Cycles', type=int, required=True, dest="Cycles",help="The number of cycles in the 4D MRI images")

	#Name of the output folder
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=True, dest="OutputFolder",help="The foldername to put all the results in.")

	args=parser.parse_args()
	ImageAnalysis4DMRIAngiography(args).Main()
	

