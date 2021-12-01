#This file will create a velocity vtu file
#from phase and magnitude data from 4D MRI.
import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import numpy as np
import vtk
import argparse
from PrincipleComponentAnalysis import PRINCIPLE_COMPONENT_ANALYSIS

class ImageAnalysisDilateMyocardium():
	def __init__(self,Args):
		self.Args=Args
		
			

	def main(self):
		#Read the vtu file
		print ("--- Reading %s"%self.Args.InputFileName)
		Volume=self.READ_VTU(self.Args.InputFileName)
	
		print ("--- Computing Princple Component Analysis and writing file")
		Centroid,Norm1,Norm2,Apex,Size=PRINCIPLE_COMPONENT_ANALYSIS(Volume).main()

		print ("--- Reading all of the centerlines coordinates")
		Centerline_Coords=self.READ_CENTERLINES_COORDS(self.Args.CenterlinesFolder)
		
		print ("--- Find point closest to apex on CL (i.e. LAD groove)")
		P_lad_groove=self.CLOSEST_POINT(Apex,Centerline_Coords)
		
		print ("--- Shift the Apex to the Lad Groove")
		Volume_Shifted=self.SHIFT_VOLUME(Volume,P_lad_groove,Apex)	

		#print ("--- Dilate the Coordinates by %d"%self.Args.Value)	
		Volume_Dilated=self.DilateVolume(Volume,Centroid,Norm1,Norm2,Apex,Size)
		
		print ("--- Write the Dilate vtu file:%s"%self.Args.OutputFileName)
		self.WRITE_VTU(self.Args.OutputFileName,Volume_Shifted)
	
	#Dilate the Myocardium
	def DilateVolume(self,Volume,Centroid,Norm1,Norm2,Apex,Size):

		#New size of the myocardium
		Size_new=(1+self.Args.Value/100.)*Size
		#Create line from centroid in both directions
		LineToApex,LineToBase=self.CreateLine(Centroid,Apex,Size)
		#Get the number of points
		Npts=Volume.GetNumberOfPoints()
		#Loop over all the points and dilate based on distance from line
		#print (dir(Volume.GetPoints().SetPoint(0,[1,2,3])))
		#exit(1)
		for i in range(0,Npts):
			coord_=Volume.GetPoint(i)
			dist_P_to_line_=np.sqrt(vtk.vtkLine.DistanceToLine(coord_,Centroid,Apex))
			dist_P_to_Apex_=np.power( np.power(coord_[0]-Apex[0],2) + np.power(coord_[1]-Apex[1],2) + np.power(coord_[2]-Apex[2],2),0.5)
			dist_Apex_to_ProjP_=np.power(np.power(dist_P_to_Apex_,2)-np.power(dist_P_to_line_,2),0.5)
			##Find the coordinate by moving from Apex a certain distance
			#coord_ProjP_=Apex+Norm1*dist_Apex_to_ProjP_
			
			#Find the normal from Centroid
			Norm_Centroid_to_P_=np.array([Centroid[0]-coord_[0],Centroid[1]-coord_[1],Centroid[2]-coord_[2]])
			Norm_Centroid_to_P_=Norm_Centroid_to_P_/np.linalg.norm(Norm_Centroid_to_P_)
			
			#Scale the coordinate so that Size increases by % value.
			scaling_=((Size-dist_Apex_to_ProjP_)/Size) * Size*(self.Args.Value/100.)
			coord_new_x_=coord_[0]-Norm_Centroid_to_P_[0]*scaling_	
			coord_new_y_=coord_[1]-Norm_Centroid_to_P_[1]*scaling_	
			coord_new_z_=coord_[2]-Norm_Centroid_to_P_[2]*scaling_	
			coord_new=np.array([coord_new_x_,coord_new_y_,coord_new_z_])
			Volume.GetPoints().SetPoint(i,coord_new)
		return Volume

	#Read all of the Centerline Coordinates
	def READ_CENTERLINES_COORDS(self,CenterlinesFolder):
		filenames=sorted(glob("%s/wall_LCA_*cl.vtp"%CenterlinesFolder))
		Coords=[]
		for filename in filenames:
			reader=vtk.vtkXMLPolyDataReader()
			reader.SetFileName(filename)
			reader.Update()
			CL_data=reader.GetOutput()
			N=CL_data.GetNumberOfPoints()
			for i in range(N):
				Coords.append(CL_data.GetPoint(i))
		return np.array(Coords)

	#Shift the volume so the apex lies on the LAD groove.
	def SHIFT_VOLUME(self,Volume,Point,Apex):
		Vector=np.array([Apex[0]-Point[0],Apex[1]-Point[1],Apex[2]-Point[2]])
		Npts=Volume.GetNumberOfPoints()
		for i in range(0,Npts):
			coord_=Volume.GetPoint(i)
			coord_=coord_-Vector
			Volume.GetPoints().SetPoint(i,coord_)
		return Volume
	
	#Create a line from apex and centroid of the myocardium
	def CreateLine(self,Centroid,Apex,Size):
		line0=np.array([Centroid[0]-Apex[0],Centroid[1]-Apex[1],Centroid[2]-Apex[2]])
		line1=-1*line0	
		line0=(line0/np.linalg.norm(line0))*(Size)
		line1=(line1/np.linalg.norm(line1))*(Size)
		return line0,line1

	def CLOSEST_POINT(self,Point, Array):
		dist_2 = np.sum((Array - Point)**2, axis=1)
		return Array[np.argmin(dist_2)]

        #Read the VTU file
	def READ_VTU(self,filename):
		Volume=vtk.vtkXMLUnstructuredGridReader()
		Volume.SetFileName(filename)
		Volume.Update()
		return Volume.GetOutput()

	def WRITE_VTU(self,filename,data):
		Volume=vtk.vtkXMLUnstructuredGridWriter()
		Volume.SetFileName(filename)
		Volume.SetInputData(data)
		Volume.Update()


if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will dilate the myocardium by a specified vaue in cm")

	parser.add_argument('-ifile', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the myocardial volume")

	parser.add_argument('-value', '--Value', type=float, required=True, dest="Value",help="The percentage of dilation as a function of Myocardial Size. Max at apex and 0 at base.")
	
	parser.add_argument('-scaling', '--Scaling', type=str, required=False, dest="Scaling",help="The method to scale from apex to base: Linear, Parabolic or Exponential")
	
        #Input folder that contains the coronary centerlines
	parser.add_argument('-centerlinesfolder', '--CenterlinesFolder', type=str, required=True, dest="CenterlinesFolder",help="Folder that contains the centerline files")

        #Output Filename 
	parser.add_argument('-ofile', '--OutputFileName', type=str, required=True, dest="OutputFileName",help="The name of the file to store the dialate LV")

        
	args=parser.parse_args()
	ImageAnalysisDilateMyocardium(args).main()
