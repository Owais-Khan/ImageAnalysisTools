#This script will compute the volume of the myocardium
#base on a given threshold. It will also check for noise
#in the image to ensure that neighbouring nodes have a positve
#value

import vtk
import numpy as np
import os
import argparse

class ImageAnalysisMyocardiumVolume():
	def __init__(self,Args):
		self.Args=Args


	def main(self):
		

        
	def ReadingDicomFolder(self,inputfoldername,outfilename):
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
	
