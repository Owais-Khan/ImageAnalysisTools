import vtk
import numpy as np
from vmtk import vmtkscripts
import os
import sys
from glob import glob
import argparse

class CoronaryCenterlines():
	def __init__(self,Args):
		self.Args=Args
		#Create a centerlines folder inside the user-defined Output Folder
		print ("Creating the output directory: %s/Centerlines/"%self.Args.OutputFolderName)
		os.system("mkdir %s/Centerlines/"%self.Args.OutputFolderName)
		#Find all of the centerlines inside the mesh-surface folder.
		filenamesLCA=sorted(glob("%s/mesh-surfaces/wall_LCA*.vtp"%self.Args.InputFolderName))
		filenamesRCA=sorted(glob("%s/mesh-surfaces/wall_RCA*.vtp"%self.Args.InputFolderName))
		self.filenames=filenamesLCA+filenamesRCA
		print ("Found %d LCA files inside the mesh-surface folder"%len(filenamesLCA))
		print ("Found %d RCA files inside the mesh-surface folder"%len(filenamesRCA))

	def main(self):
		for filename in self.filenames:
			print ("----Computing Centerlines for: %s"%filename)
			#Read the surface
			surface=self.VTP_Reader(filename)
			#Compute Centerlines
			centerline=vmtkscripts.vmtkCenterlines()
			centerline.Surface=surface
			centerline.SeedSelectorName="openprofiles"
			centerline.Interactive=1
			centerline.Resampling=1
			centerline.ResamplingStepLength=0.06
			centerline.Execute()
			#Write the centerline file
			surfaceWriter=vmtkscripts.vmtkSurfaceWriter()
			surfaceWriter.Surface=centerline.Centerlines
			surfaceWriter.OutputFileName="%s/Centerlines/%s"%(self.Args.OutputFolderName,filename.split("/")[-1].replace(".vtp","_cl.vtp"))
			surfaceWriter.Execute()
	def VTP_writer(self,centerline,filename):
		filename_new="%s/Centerlines/%s"%(self.Args.OutputFolderName,filename.split("/")[-1].replace(".vtp","_cl.vtp"))
		print ("----Writing: %s",filename_new)
		surfaceWriter=vtk.vtkXMLPolyDataWriter()
		surfaceWriter.SetFileName(filename_new)
		surfaceWriter.SetInputData(centerline)
		surfaceWriter.Update()

	
	def VTP_Reader(self,filename):
		surfaceReader=vtk.vtkXMLPolyDataReader()
		surfaceReader.SetFileName(filename)
		surfaceReader.Update()
		surface=surfaceReader.GetOutput()
		return surface

if __name__=='__main__':
        #Description
	parser = argparse.ArgumentParser(description="This script will compute the centerlines of coronary arteries using the 'mesh-complete' that simvascular outputs")

        #Input filename of the perfusion map
	parser.add_argument('-ifolder', '--InputFolderName', type=str, required=True, dest="InputFolderName",help="The path to the mesh-complete folder")

        #Resolution of the Averaging Procedure
	parser.add_argument('-ofolder', '--OutputFolderName', type=str, required=True, dest="OutputFolderName",help="The folder path for the output folder that will contain centerline files")


	args=parser.parse_args()
	CoronaryCenterlines(args).main()
	
