import vtk
import numpy as np
from utilities import *
from glob import glob
import argparse

class ImageAnalysisPickRandomPoints():
	def __init__(self,Args):
		self.Args = Args
		self.Args.InputFiles=sorted(glob(self.Args.InputFolder+"/*.vtu"))
		
		#Initialze arrays to store the coordinates
		self.X_sen=np.zeros(self.Args.NumberOfPoints)
		self.Y_sen=np.zeros(self.Args.NumberOfPoints)
		self.Z_sen=np.zeros(self.Args.NumberOfPoints)
		self.Vx_sen=np.zeros(self.Args.NumberOfPoints)
		self.Vy_sen=np.zeros(self.Args.NumberOfPoints)
		self.Vz_sen=np.zeros(self.Args.NumberOfPoints)

	def Main(self):
		#Initialize arrays to store the coordinates	
		counter=0
		for FileName_ in self.Args.InputFiles:
			print ("Looping over:%s "%FileName_)
			VelocityData_=ReadVTUFile(FileName_)
			if counter==0:
				X,Y,Z,Vx,Vy,Vz=self.GetNonZeroVelocityCoords(VelocityData_)
				
			#Pick random coordinates
			 
	def GetNonZeroVelocityCoords(self,VelocityData):
		Npts=VelocityData.GetNumberOfPoints()
		X_nonzero=[]; Vel_X=[]
		Y_nonzero=[]; Vel_Y=[]
		Z_nonzero=[]; Vel_Z=[]
		counter1=0
		tol=1e-8
		for i in range(0,Npts*3,3):
			Vel_x_=VelocityData.GetPointData().GetArray("velocity").GetValue(i)
			Vel_y_=VelocityData.GetPointData().GetArray("velocity").GetValue(i+1)
			Vel_z_=VelocityData.GetPointData().GetArray("velocity").GetValue(i+2)
			if Vel_x_<tol and Vel_y_<tol and Vel_z_<tol:
				continue
			else:
				X_nonzero.append(VelocityData.GetPoint(counter1)[0])
				Y_nonzero.append(VelocityData.GetPoint(counter1)[1])
				Z_nonzero.append(VelocityData.GetPoint(counter1)[2])
				Vel_X.append(Vel_x_)
				Vel_Y.append(Vel_y_)
				Vel_Z.append(Vel_z_)
			counter1+=1

		return X_nonzero,Y_nonzero,Z_nonzero, Vel_X,Vel_Y,Vel_Z		


if __name__ == "__main__":
    	#Description
	parser = argparse.ArgumentParser(description = "This script will take a folder that contains volume files and pick random points inside the domain. It looks for velocity vector field and ensures only coordinates with non-zero velocity, away from the wall, are picked.")

    	#InputFiles
	parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = "InputFolder", help = "The folder that contains the velocity volumetric files.")
	
	parser.add_argument('-NumberOfPoints', '--NumberOfPoints', type = int, required = True, dest = "NumberOfPoints", help = "The number of points to extract from the volumetric files")

	args = parser.parse_args()
	ImageAnalysisPickRandomPoints(args).Main()



