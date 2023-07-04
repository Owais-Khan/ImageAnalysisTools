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
		
		if self.Args.OutputFolderName is None:
			if self.Args.InputFolderName[-1]=="/": 
				self.Args.OutputFolderName=self.Args.InputFolderName[:-1]+"_Centerlines"
			else:
				self.Args.OutputFolderName=self.Args.InputFolderName+"_Centerlines"
				
		

		#Create a centerlines folder inside the user-defined Output Folder
		print ("Creating the output directory: %s"%self.Args.OutputFolderName)
		os.system("mkdir %s/"%self.Args.OutputFolderName)
		#Find all of the centerlines inside the mesh-surface folder.
		filenamesLCA=sorted(glob("%s/L_*.vtp"%self.Args.InputFolderName))
		filenamesRCA=sorted(glob("%s/R_*.vtp"%self.Args.InputFolderName))
		self.filenames=filenamesLCA+filenamesRCA
		print ("Found %d LCA files inside the input folder"%len(filenamesLCA))
		print ("Found %d RCA files inside the input folder"%len(filenamesRCA))

	def main(self):
		#Automatically estimate the centerline resample rate
		print ("--- Reading L_LAD.vtp surface to estimate CL resolution")
		surface=self.VTP_Reader(self.Args.InputFolderName+"/L_LAD.vtp")
		#Get minimum and max data ranges
		BBox=surface.GetBounds()
		Length=((BBox[1]-BBox[0])**2+(BBox[3]-BBox[2])**2+(BBox[5]-BBox[4])**2)**0.5
		if Length<15: CL_Res=0.06
		else: CL_Res=0.6
		print ("--- Bounding Box of LAD is: %.05f"%Length)
		print ("--- Centerline Resolution: %.05f"%CL_Res)

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
			centerline.ResamplingStepLength=CL_Res
			centerline.Execute()
			#Write the centerline file
			surfaceWriter=vmtkscripts.vmtkSurfaceWriter()
			surfaceWriter.Surface=centerline.Centerlines
			surfaceWriter.OutputFileName="%s/%s"%(self.Args.OutputFolderName,filename.split("/")[-1].replace(".vtp","_cl.vtp"))
			surfaceWriter.Execute()
	def VTP_writer(self,centerline,filename):
		filename_new="%s/%s"%(self.Args.OutputFolderName,filename.split("/")[-1].replace(".vtp","_cl.vtp"))
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
	parser = argparse.ArgumentParser(description="This script will compute the centerlines of coronary arteries using the vessel labels you generate from SimVascular. R=Right, L=Left")

        #Input filename of the perfusion map
	parser.add_argument('-InputFolderName', '--InputFolderName', type=str, required=True, dest="InputFolderName",help="The path containing the vessel labels .")

        #Resolution of the Averaging Procedure
	parser.add_argument('-OutputFolderName', '--OutputFolderName', type=str, required=False, dest="OutputFolderName",help="The folder path for the output folder that will contain centerline files. If not, select 'Results'")


	args=parser.parse_args()
	CoronaryCenterlines(args).main()
	
