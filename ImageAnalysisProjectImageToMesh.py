import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import os
from glob import glob
import argparse
from utilities import *
class ImageAnalysisProjectImageToMesh():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFileName is None:
			self.Args.OutputFileName=self.Args.InputFileName2[0:-4]+"_"+self.Args.InputFileName1.split("/")[-1][0:-4]+self.Args.InputFileName2[-4:]
			print ("--- OutputFileName not provided...")
			print ("--- Using %s"%self.Args.OutputFileName)
			
	def Main(self):
		#Read the source data
		print ("--- Loading the source data: %s"%self.Args.InputFileName1)
		if self.Args.InputFileName1[-4:]==".vti":
			SourceData=ReadVTIFile(self.Args.InputFileName1)
		elif self.Args.InputFileName1[-4:]==".vtu":
			SourceData=ReadVTUFile(self.Args.InputFileName1)     
		elif self.Args.InputFileName1[-4:]==".vtp":
			SourceData=ReadVTPFile(self.Args.InputFileName1) 
		else:
			print ("Input file format not detected. Exiting...")
			exit(1)

		print ("--- Loading the target volume/surface mesh: %s"%self.Args.InputFileName2)
                
		if self.Args.InputFileName2[-4:]==".vtu":
			TargetData=ReadVTUFile(self.Args.InputFileName2)
		elif self.Args.InputFileName2[-4:]==".vtp":
			TargetData=ReadVTPFile(self.Args.InputFileName2)
		else:
			print ("Target mesh format not detected. Exiting...")
			exit(1)	


		print ("--- Create a probe filter to interpolate source to target")
		ProbeFilter=vtk.vtkProbeFilter()
		ProbeFilter.SetInputData(TargetData)
		ProbeFilter.SetSourceData(SourceData)
		ProbeFilter.Update()
		ProbeOutput=ProbeFilter.GetOutput()
	
		#Loop over the surface nodes and write a velocity of zero
		print ("--- Setting Velocity on the Surface to Zero")
		if self.Args.InputFileName3 is not None:
			print ("--- Reading Wall Surface Mesh where Zero Velocity is Assigned")
			Surface=ReadVTPFile(self.Args.InputFileName3)
		else:
			print ("--- No Wall Mesh Found. Assign zero velocity on walls and caps")
			Surface=ExtractSurface(TargetData)


		SurfaceCoords=vtk_to_numpy(Surface.GetPoints().GetData())
		counter=0
		for i in range(ProbeOutput.GetNumberOfPoints()):
			if TargetData.GetPoint(i) in SurfaceCoords:
				ProbeOutput.GetPointData().GetArray(self.Args.ArrayName).SetValue(i*3,0)
				ProbeOutput.GetPointData().GetArray(self.Args.ArrayName).SetValue(i*3+1,0)
				ProbeOutput.GetPointData().GetArray(self.Args.ArrayName).SetValue(i*3+2,0)
		
	
	
		print ("--- Write the output file in the same format as target data")
		WriteVTUFile(self.Args.OutputFileName,ProbeOutput)		
		

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will extract the image intensities and interpolate them onto a volumetric mesh.")

	parser.add_argument('-InputFileName1', '--InputFileName1', type=str, required=True, dest="InputFileName1",help="File name of the source data")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, default="Velocity",dest="ArrayName",help="The array name that contains the image intensitites.")
        
	#Input filename of the coronary segmented surface.
	parser.add_argument('-InputFileName2', '--InputFileName2', type=str, required=True, dest="InputFileName2",help="File name of the input mesh/surface on which to project from the source data")
	
	parser.add_argument('-InputFileName3', '--InputFileName3', type=str, required=False, dest="InputFileName3",help="File name of the input wall mesh where velocity of zero should be assigned")
        
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="File name in which to store the image intensitites.")
        
	args=parser.parse_args()
       
	ImageAnalysisProjectImageToMesh(args).Main()

