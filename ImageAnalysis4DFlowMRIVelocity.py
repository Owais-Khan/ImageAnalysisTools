import numpy as np
import vtk
import argparse
from glob import glob
from utilities import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class ImageAnalysis4DFlowMRIVelocity():
	def __init__(self,Args):
		self.Args=Args
		
		if self.Args.OutputFolder is None:
			OutputFolder=self.Args.InputFolder.replace(self.Args.InputFolder.split("/")[-1],"VelocityData")
			print ("---The OutputFolder is: %s"%OutputFolder)	


	def Main(self):
		#Read the Mesh 
		Mesh=ReadVTUFile(self.Args.InputFile)
		
		#Read all of the Phase Files
		Phase1Files=[]
		Phase2Files=[]
		Phase3Files=[]
		AllFiles=sorted(glob(self.Args.InputFolder+"/Phase*.vti"))
		for FileName_ in AllFiles:
			if FileName_.find("Phase1")>=0:
				Phase1Files.append(FileName_)
			if FileName_.find("Phase2")>=0:
				Phase2Files.append(FileName_)
			if FileName_.find("Phase3")>=0:
				Phase3Files.append(FileName_)

		print ("---The number of Phase 1 Files are: %d"%len(Phase1Files))
		print ("---The number of Phase 2 Files are: %d"%len(Phase2Files))
		print ("---The number of Phase 3 Files are: %d"%len(Phase3Files))
		
		#Loop over all of the phase files
		for i in range(0,len(Phase1Files)):
			print ("--- Looping over Files %s of %s"%(i,len(Phase1Files)))
			
			Phase1FileData_=ReadVTIFile(Phase1Files[i])
			print ("------ Read: %s"%Phase1Files[i])	
			Phase2FileData_=ReadVTIFile(Phase2Files[i])	
			print ("------ Read: %s"%Phase2Files[i])	
			Phase3FileData_=ReadVTIFile(Phase3Files[i])	
			print ("------ Read: %s"%Phase3Files[i])	
			
			#Project Image Intensities onto Volume Mesh
			Phase1Intensities_=self.ProbeFilter(Mesh,Phase1FileData_)
			Phase2Intensities_=self.ProbeFilter(Mesh,Phase2FileData_)
			Phase3Intensities_=self.ProbeFilter(Mesh,Phase3FileData_)
	
			if i==0:
				Npts=Phase1Intensities_.GetNumberOfPoints()
				Velocity=np.zeros(shape=(Npts,3))
			for j in range(0,Npts):
				Velocity[j,0]=Phase1Intensities_.GetPointData().GetArray(self.Args.ImageArrayName).GetValue(j)/self.Args.ImageArrayScaling*self.Args.Venc
				Velocity[j,1]=Phase2Intensities_.GetPointData().GetArray(self.Args.ImageArrayName).GetValue(j)/self.Args.ImageArrayScaling*self.Args.Venc
				Velocity[j,2]=Phase3Intensities_.GetPointData().GetArray(self.Args.ImageArrayName).GetValue(j)/self.Args.ImageArrayScaling*self.Args.Venc

			numpy_to_vtk	


	def ProbeFilter(self,TargetData,SourceData):
		ProbeFilter=vtk.vtkProbeFilter()
		ProbeFilter.SetInputData(TargetData)
		ProbeFilter.SetSourceData(SourceData)
		ProbeFilter.Update()
		ProbeOutput=ProbeFilter.GetOutput()
		return ProbeOutput
		
		

			
			
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take a volumetric mesh and projct 4D Flow MRI velocities onto the mesh.")

	parser.add_argument('-InputFile', '--InputFile', required=True, dest="InputFile",help="A .vtu mesh file.")
	
	parser.add_argument('-InputFolder', '--InputFolder', required=True, dest="InputFolder",help="A folder that contains files with Phase 1, 2 and 3 tags.")
	
	parser.add_argument('-Venc', '--Venc', required=False, default=150, dest="Venc",help="This velocity encoding parameter from the MRI scan.")
	
	parser.add_argument('-ScalingFactor', '--ScalingFactor', required=False, default=1, dest="ScalingFactor",help="Used to scale the mesh from mm to cm if needed.")
	
	parser.add_argument('-ImageArrayName', '--ImageArrayName', required=False, default="DICOMImage", dest="ImageArrayName",help="The array name of the pixel intensities stored in the Phase 1, 2 and 3 files.")
	
	parser.add_argument('-ImageArrayScaling', '--ImageArrayScaling', required=False, default=4096, dest="ImageArrayScaling",help="The scaling factor used to scale the image array.")
	
	parser.add_argument('-OutputFolder', '--OutputFolder', required=False, dest="OutputFolder",help="The folder name where all the output velocities will be stored.")



	args=parser.parse_args()
	ImageAnalysis4DFlowMRIVelocity(args).Main()
	

