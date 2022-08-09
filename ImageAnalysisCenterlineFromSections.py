#THis code was writin by Sensi Vivian on July 28
# Curvature Projection
# July 28, 2022

import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import vtk
import numpy as np
import argparse
from utilities import *

class CenterlineFromSections():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		if self.Args.OutputFolder == None:
			self.Args.OutputFolder="Results"
			os.system("mkdir %s"%self.Args.OutputFolder)

		outfile=open(self.Args.OutputFolder+"/CenterlineFromSections.txt",'w')
		outfile.write('SectionID,Ox,Oy,Oz\n')

		CenterlineSections={}
		CenterlineGeometry={}
		CLSection_={}
		file = self.Args.InputSurface.split('.')[0]
		#Extract the centerlines
		"""os.system("vmtkcenterlines -ifile %s.vtp -ofile %s_CL.vtp -endpoints 1 -resampling 1 -resamplingstep 0.05"%(file,file))

		#Convert surface to mesh 
		os.system("vmtksurfacetomesh -ifile %s.vtp -ofile %s.vtu"%(file,file))

		#Compute centerline sections
		os.system("vmtkcenterlinemeshsections -centerlinesfile %s_CL.vtp -ifile %s.vtu -ofile %s_CLsections.vtp"%(file,file,file))"""

		CenterlineSections=ReadVTPFile(file+"_CLsections.vtp")
		 
		Nstart,Nend=CenterlineSections.GetPointData().GetArray("SectionIds").GetRange()
		
		#Convert to an integer
		Nstart=int(Nstart)
		Nend=int(Nend)

	
		#Get the distance and centroid
		Distance=np.zeros(Nend)
		Centroid=[]

		for i in range (0,Nend):
			#Get the Section of the Slice
			section_=ThresholdInBetween(CenterlineSections,"SectionIds",i,i)
				
			#Extract a surface from the section
			surface_=vtk.vtkDataSetSurfaceFilter()
			surface_.SetInputData(section_)
			surface_.Update()
			
			#Get the centroid
			Centroid.append(list(GetCentroid(surface_.GetOutput())))
				
			if i>0:
				Distance[i]=np.sqrt( (Centroid[i][0]-Centroid[i-1][0])**2 + (Centroid[i][1]-Centroid[i-1][1])**2 + (Centroid[i][2]-Centroid[i-1][2])**2)	
			outfile.write("%d %0.5f %0.5f %.05f %.05f\n"%(i,Centroid[i][0],Centroid[i][1],Centroid[i][2],Distance[i]))
				
		outfile.close()	

		#Create a centerline using the centroid points
		PolyLine=ConvertPointsToLine(Centroid)	

		#Write the CL as surface file
		WriteVTPFile("%s/CenterlineFromMeshSections.vtp"%self.Args.OutputFolder,PolyLine)	
	
	
if __name__=="__main__":
	#Description
	parser = argparse.ArgumentParser(description="This script will compute a centerline based on the cross-section of the model")
	parser.add_argument('-InputSurface','--InputSurface', type=str, required=True, dest="InputSurface", help="The surface for which we need the centerline")

	parser.add_argument('-OutputFolder','--OutputFolder',type=str, required=False, dest="OutputFolder",help="An output folder to store the reigstered surfaces and output results")
	args=parser.parse_args()
	CenterlineFromSections(args).Main()
