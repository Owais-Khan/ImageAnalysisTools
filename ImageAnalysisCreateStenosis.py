import vtk
import numpy as np
from vmtk import vmtkscripts
import os
import sys
from glob import glob
import argparse
from utilities import *

class ImageAnalysisCreateStenosis():
	def __init__(self,Args):
		self.Args=Args
	
	def Main(self):
		#Read the Input Surface File
		print ("--- Reading the Input Surface File: %s"%self.Args.InputFileName)
		InputSurface=ReadVTPFile(self.Args.InputFileName)
		
		#Read the Input Centerline File
		print ("--- Reading the Input Centerlines File: %s"%self.Args.InputCenterlinesFileName)
		InputCenterlines=ReadVTPFile(self.Args.InputCenterlinesFileName)
		
		#Get the length and MSIR array from centerlines
		print ("--- Computing Centerline Parameters (e.g., length, radii)")
		Length,MISR,CLxyz=self.GetCenterlineParameters(InputCenterlines)
		
		#Get the location of the stenosis
		print ("--- Generating the stenosis function")
		LocStenosis=Length[-1]*(self.Args.StenosisLocation/100.) #stenosis location
		MISR_avg=np.average(MISR)
		Dia=2*MISR_avg #Average diameter of the aorta
		s0=-0.5*self.Args.StenosisDiameter+0.5
		LengthStenosis=2*Dia
		DiameterFactor=np.zeros(InputCenterlines.GetNumberOfPoints())
		x_=0	
		#Stenosis equation from Varghese 2007 J. Fluid Mechanics
		for i in range(0,InputCenterlines.GetNumberOfPoints()):	
			if (Length[i]>=(LocStenosis-Dia)) and (Length[i]<=(LocStenosis+Dia)):
				x_+=Length[i]-Length[i-1]
				S_=0.5*Dia*(1-s0*(1+np.cos(2.*np.pi*(Length[i]-LocStenosis)/LengthStenosis)))
				DiameterFactor[i]=S_/Dia*2
			else:
				DiameterFactor[i]=1.0

		StenosisMinDia=min(DiameterFactor)*100
		print ("------ Stenosis Throat has a %.01f %% of Normal Diameter"%(StenosisMinDia))

		#Loop over all of the surface nodes, find closest centerline point, and reduce diameter
		print ("--- Updating the Input Surface to Include the Stenosis")
		for i in range(0,InputSurface.GetNumberOfPoints()):
			ClosestCLPoint,ClosestCLPointArg=ClosestPoint(np.array(InputSurface.GetPoint(i)),np.array(CLxyz))
			if DiameterFactor[ClosestCLPointArg]==1.0:
				continue
			else:
				#Get Normal array between SurfacePoint and CLPoint
				CLPoint_=CLxyz[ClosestCLPointArg] #Closest Centerline Point
				SurfPoint_=InputSurface.GetPoint(i) #Current Surface Point
				NormX_=(CLPoint_[0]-SurfPoint_[0])*(1-DiameterFactor[ClosestCLPointArg])
				NormY_=(CLPoint_[1]-SurfPoint_[1])*(1-DiameterFactor[ClosestCLPointArg])
				NormZ_=(CLPoint_[2]-SurfPoint_[2])*(1-DiameterFactor[ClosestCLPointArg])

				#Get New Surface Point
				SurfPointNew_=np.array([SurfPoint_[0]+NormX_,SurfPoint_[1]+NormY_,SurfPoint_[2]+NormZ_])
				InputSurface.GetPoints().SetPoint(i,SurfPointNew_)

			
		#Write the output surface
		self.Args.OutputFileName=self.Args.InputFileName.replace(".vtp","_stenosis%d.vtp"%(self.Args.StenosisDiameter*100))
		print ("--- Writing the Output File: %s"%self.Args.OutputFileName)
		WriteVTPFile(self.Args.OutputFileName,InputSurface)
			

	def GetCenterlineParameters(self,Centerlines):
		Length=np.zeros(Centerlines.GetNumberOfPoints())
		MISR=np.zeros(Centerlines.GetNumberOfPoints())
		CenterlineCoordinates=[Centerlines.GetPoint(0)]
		StenosisTag=[]
		for i in range(1,Centerlines.GetNumberOfPoints()):
			CenterlineCoordinates.append(Centerlines.GetPoint(i))
			xyz0=Centerlines.GetPoint(i-1)
			xyz1=Centerlines.GetPoint(i)
			MISR[i]=Centerlines.GetPointData().GetArray("MaximumInscribedSphereRadius").GetValue(i)
			distance=((xyz0[0]-xyz1[0])**2+(xyz0[1]-xyz1[1])**2+(xyz0[2]-xyz1[2])**2)**0.5
			Length[i]=Length[i-1]+distance

	
		return Length,MISR,CenterlineCoordinates 


if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will artifically create a stenosis at a prescribe length along the centerline.")

        #Input filename of the perfusion map
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The file name of the surface that contains the geoemtry.")

	parser.add_argument('-InputCenterlinesFileName', '--InputCenterlineFileName', type=str, required=True, dest="InputCenterlinesFileName",help="The file name of the centerlines that contains the centerline of the surface geometry. No branch centerlines should be included, only the centerline of the artery where a stenosis is needed.")
	
	parser.add_argument('-StenosisLocation', '--StenosisLocation', type=float, required=True, dest="StenosisLocation",help="The location as a percentage of the length of the centerline. For example, 50 would mean create a stenosis at 50% of the length of the centerline.")
	
	parser.add_argument('-StenosisLength', '--StenosisLength', type=float, required=False, default=2, dest="StenosisLength",help="The length of the stenosis in terms of maximum inscripted sphere radius (MISR). For example, StenosisLength of 1 would mean that stenosis is 0.5*MISR proximal and 0.5*MISR distal to the stenosis throat. Default is 2.0")
	
	parser.add_argument('-StenosisDiameter', '--StenosisDiameter', type=float, required=False, default=0.5, dest="StenosisDiameter", choices=np.round(np.linspace(0,1.,101),decimals=2), help="The minimum diameter of the stenosis as a fraction of the global diameter. Default is 0.5")

        #Resolution of the Averaging Procedure
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The filename of the output surface file with the stenosis included.")

	args=parser.parse_args()
	ImageAnalysisCreateStenosis(args).Main()
