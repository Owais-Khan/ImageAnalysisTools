# Centerline Deviation Script
# Resamples and computes the difference between two clipped and registered centerlines
# August 19, 2022

import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import vtk
import numpy as np
import argparse
from utilities import *

class CenterlineDeviation():
	def __init__(self,Args):
                self.Args=Args

	def Main(self):
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder="CenterlineDeviation_Results"
			os.system("mkdir %s"%self.Args.OutputFolder)
		
		Centerlines = [self.Args.BaseCenterline,self.Args.InputCenterline]
		Length=[]
		for filename in Centerlines:
			Dist_=0
			file = ReadVTPFile(filename)

			#Compute length of lines
			for i in range (1,file.GetNumberOfPoints()):
				Dist_+=np.sqrt((file.GetPoint(i)[0]-file.GetPoint(i-1)[0])**2+(file.GetPoint(i)[1]-file.GetPoint(i-1)[1])**2+(file.GetPoint(i)[2]-file.GetPoint(i-1)[2])**2)
			Length.append(Dist_)	
			print ("--- The Length of %s is: %.04f"%(filename,Dist_))
		
		#Resample lines to have equal amount of points
		resampleLength=Length[0]/self.Args.NPoints  
		os.system("vmtkcenterlineresampling -ifile %s -length %f -ofile %s"%(Centerlines[0],resampleLength,self.Args.OutputFolder+'/'+Centerlines[0].split('/')[-1]))

		#Read Resampled Centerlines
		BaseCL_=ReadVTPFile(self.Args.OutputFolder+'/'+Centerlines[0].split('/')[-1])
		InputCL_=ReadVTPFile(Centerlines[1])
		BaseCLPts_=np.array([BaseCL_.GetPoint(i) for i in range(BaseCL_.GetNumberOfPoints())])
		InputCLPts_=np.array([InputCL_.GetPoint(i) for i in range(InputCL_.GetNumberOfPoints())])
	
		Mag_=np.zeros(BaseCL_.GetNumberOfPoints())
		PolyData=vtk.vtkAppendPolyData()

		outfile = open ("CenterlineDeviation_Results.txt",'w')
		outfile.write('PointID,Distance,Deviation_Magnitude\n')
		
		#Compute difference between centerlines	
		Length_=0
		LengthIncrement_=0
		for i in range (BaseCL_.GetNumberOfPoints()):
			ClosestPoint_,Dist_=ClosestPoint(BaseCLPts_[i],InputCLPts_)
			points = [BaseCLPts_[i], ClosestPoint_]
		
			if i>0:
				LengthIncrement_=np.sqrt(sum((BaseCLPts_[i]-BaseCLPts_[i-1])**2))
				Length_+=LengthIncrement_
		
			if Length_<Length[1]:
				Mag_[i]=np.sqrt((points[0][0]-points[1][0])**2+(points[0][1]-points[1][1])**2+(points[0][2]-points[1][2])**2)
					
				PolyData.AddInputData(ConvertPointsToLine(points))
				PolyData.Update()
			else:
				print ("The Max Length of Input Centerline Reached!!!")
				print ("Stopping...")
				break
		
			outfile.write("%d,%.05f,%.05f\n"%(i,Length_,Mag_[i]))

		#Add array and write file	
		WriteVTPFile("%s/CenterlineDeviation.vtp"%self.Args.OutputFolder,PolyData.GetOutput())
		CenterlineData=ReadVTPFile("%s/CenterlineDeviation.vtp"%self.Args.OutputFolder)

		MagVTK=numpy_to_vtk(Mag_)
		MagVTK.SetName("Magnitude")
		CenterlineData.GetCellData().AddArray(MagVTK)

		WriteVTPFile("%s/CenterlineDeviation.vtp"%self.Args.OutputFolder,CenterlineData)

if __name__=="__main__":
	#Description
	parser = argparse.ArgumentParser(description="This script will compute the deviation between two registered centerlines")
	parser.add_argument('-InputCenterline','--InputCenterline', type=str, required=True, dest="InputCenterline", help="The centerline for which we will be comparing to the base centerline")
	parser.add_argument('-BaseCenterline','--BaseCenterline', type=str, required=True, dest="BaseCenterline",help="The reference centerline")
	parser.add_argument('-NPoints','--NPoints', type=int, required=False, default=200, dest="NPoints", help="Number of points used to resample and compute the centerline deviation")

	parser.add_argument('-OutputFolder','--OutputFolder',type=str, required=False, dest="OutputFolder",help="An output folder to store the registered surfaces and output results")

	args=parser.parse_args()
	CenterlineDeviation(args).Main()
