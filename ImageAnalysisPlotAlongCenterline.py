import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import numpy as np
import vtk
import argparse
from utilities import *

class ImageAnalysisPlotAlongCenterline():
	def __init__(self,Args):
		self.Args=Args
		

	def Main(self):
		#Compute the centerlines
		CL_FileName=self.Args.OutputFolder+"/"+self.Args.InputSurface.replace(".vtp","_cl.vtp").split("/")[-1]
		
		os.system("vmtkcenterlines -ifile %s -ofile %s -endpoints 1 -resampling 1 -resamplingstep 0.05"%(self.Args.InputSurface,CL_FileName))
		
		
		#Compute centerline sections
		OutputFileName=self.Args.OutputFolder+"/"+self.Args.InputVolume.split("/")[-1].replace(".vtu","_sections")
		os.system("vmtkcenterlinemeshsections -centerlines %s -ifile %s -sectionpointsfile %s.vtp -ofile %s.dat"%(CL_FileName,self.Args.InputVolume,OutputFileName,OutputFileName))
		
		self.InputSurface=ReadVTPFile(OutputFileName+".vtp")


		OutputFileNameData = open(OutputFileName+".txt",'w')
        	outfile.write('SectionID,Length,FlowRate,Velocity_Max,Velocity_Average\n')

		#Get the number of section ids
		Nstart,Nend=InputSurface.GetPointData().GetArray("SectionIds").GetRange()
		Nstart=int(Nstart)
		Nend=int(Nend)


		for i in range(Nstart,Nend):
        		section_ = ThresholdInBetween(InputSurface,"SectionIds",i,i)
        		Vmin_,Vmax_ =section_.GetPointData().GetArray("velocity_mag_average").GetValueRange()
        		
			#Compute the area
			
			

			#Vmax.extend([Vmax_])
        		Vsection_ = []
        		Nstart_,Nend_ =section_.GetPointData().GetArray("velocity_mag_average").GetRange()
        		for j in range(int(Nstart_),int(Nend_)):
                		Vsection_.extend([section_.GetPointData().GetArray("velocity_mag_average").GetValue(j)])
        		Vavg_ = sum(Vsection_)/len(Vsection_)
        		Vmax_ = max(Vsection_)
        		#Vavg.extend([sum(Vavg_)/len(Vavg_)])
        		with open(vel,'a') as writefile:
                		writefile.write('%s, %s, %s\n'%(i,Vmax_,Vavg_))

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute quantities along the centerline")
	parser.add_argument('-InputVolume', '--InputVolume', type=str, required=True, dest="InputVolume",help="The vtu file that contains the data in vtu format")
	
	parser.add_argument('-InputSurface', '--InputSurface', type=str, required=True, dest="InputSurface",help="The surface file along which to compute the data.")
        
	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="An output folder to store the centerlines and data")
	args=parser.parse_args()
	ImageAnalysisPlotAlongCenterline(args).Main()	
