import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *

class ImageAnalysisLabelInfarctedRegion():
	def __init__(self,Args):
		self.Args=Args
		self.Args.OutputFileName=self.Args.InputFile1.replace(".vtp","_infarctlabel.vtp")

	def Main(self):
		#Read the Myocardial Surface
		MyoSurface=ReadVTPFile(self.Args.InputFile1)

		#Read the infarcted region surface
		InfarctedSurface=ReadVTPFile(self.Args.InputFile2)

		#Get the number of points for myocardium and infarct region
		Npts_Infarct=InfarctedSurface.GetNumberOfPoints()
		
		for i in range(Npts_Infarct):
			PointId_=InfarctedSurface.GetPointData().GetArray("vtkOriginalPointIds").GetValue(i)
			MyoSurface.GetPointData().GetArray("MBF_Normalized50Q").SetValue(PointId_,self.Args.InfarctValue)
			MyoSurface.GetPointData().GetArray("MBF_Normalized75Q").SetValue(PointId_,self.Args.InfarctValue)
			MyoSurface.GetPointData().GetArray("MBF_NormalizedAVG").SetValue(PointId_,self.Args.InfarctValue)
			MyoSurface.GetPointData().GetArray("MBF_WallAveraged").SetValue(PointId_,self.Args.InfarctValue)
		#Write the output file		
		WriteVTPFile(self.Args.OutputFileName,MyoSurface)


if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will label the infarcted region (InputFile1) on the myocardium (InputFile2) with a default value of -1000. This value can be changed using the parameter InfarctValue")
        
	parser.add_argument('-InputFile1', '--InputFile1', type=str, required=True, dest="InputFile1",help="The myocardial surface for which infarct region needs to be labelled.")

	parser.add_argument('-InputFile2', '--InputFile2', type=str, required=True, dest="InputFile2",help="The infarcted surface extracted from manually cropping in Paraview using 'weights' (average number of nodes across the myocardial thickness).")
        
	parser.add_argument('-InfarctValue', '--InfarctValue', type=float, required=False, default=-999, dest="InfarctValue",help="The numerical value to assign to the infarcted region.")

        #Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to store the output file")

	args=parser.parse_args()
	ImageAnalysisLabelInfarctedRegion(args).Main()
 
