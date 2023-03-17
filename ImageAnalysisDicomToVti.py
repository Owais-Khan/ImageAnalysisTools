import numpy as np
import vtk
import argparse
from glob import glob

class ImageAnalysisDicomToVti():
	def __init__(self,Args):
		self.Args=Args
		self.Args.InputFolders=sorted(glob(self.Args.InputFolder+"/*"))
	def Main(self):
		#Loop over all of the input folders
		for InputFolder_ in self.Args.InputFolders:
			print ("--- Looping over %s"%InputFolder_)
			if InputFolder_[-1] is not "/": InputFolder_+="/"
			self.ConvertDicomDirectoryToVti(InputFolder_)
        
	def ConvertDicomDirectoryToVti(self,InputFolder):
		#Output FileName
		OutputFileName=self.Args.OutputFolder+InputFolder.split("/")[-2]
		print ("--- Writing %s"%OutputFileName)
		
		#Read the dicom files in the directory
		Data=self.ReadDicomDirectory(InputFolder)
		Spacing=np.array(Data.GetPixelSpacing())*self.Args.ScalingFactor
		Origin=np.array(Data.GetImagePositionPatient())*self.Args.ScalingFactor
		Data=Data.GetOutput()
			
		#Set the Origin and Spacing
		Data.SetSpacing(Spacing)
		Data.SetOrigin(Origin)
		
		#Write the Dicom File
		if self.Args.OutputFormat=="vtk":
			self.WriteData(Data,OutputFileName+".vtk")
		else:
			self.WriteVtiFile(Data,OutputFileName+".vti")

	def ReadVtiFile(self,FileName):
		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(FileName)
		reader.Update()
		return reader.GetOutput()

	def ReadDicomDirectory(self,FolderName):
		reader=vtk.vtkDICOMImageReader()
		reader.SetDirectoryName(FolderName)
		reader.Update()
		return reader

	def WriteData(self,Data,FileName):
		writer=vtk.vtkDataSetWriter()
		writer.SetFileName(FileName)
		writer.SetInputData(Data)
		writer.Update()

	def WriteVtiFile(self,Data,FileName):	
		#Write the file in vti format
		writer=vtk.vtkXMLImageDataWriter()
		writer.SetFileName(FileName)
		writer.SetInputData(Data)
		writer.Write()
	

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take input folder(s) and write a vti file in each folder")

        #Provide a path to the Magnitude Images
	parser.add_argument('-InputFolder', '--InputFolder', required=True, dest="InputFolder",help="A folder that contains a bunch of sub-folders each with their own image stacks")
	
	parser.add_argument('-OutputFolder', '--OutputFolder', required=True, dest="OutputFolder",help="The folder name where the script stores all of the .vti and vtk image stacks.")
	
	parser.add_argument('-OutputFormat', '--OutputFormat', required=False, dest="OutputFormat", default="vti", help="The output format for the image stack. vti by default otherwise, you may specify vtk")

	#Provide a scaling factor	
	parser.add_argument('-ScalingFactor', '--ScalingFactor', type=float, required=False, default=1.0, dest="ScalingFactor", help="Scaling factor to convert images from one unit to another")

        
	args=parser.parse_args()
	ImageAnalysisDicomToVti(args).Main()
 
