#This script will read the territories file
#and write a tecplot file with histogram plots

import vtk
import numpy as np
import scipy as sp
import os
from glob import glob
from scipy.spatial import distance as DISTANCE
from scipy import stats
import argparse
from utilities import *
import copy
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class CoronaryPerfusionPlot():
	def __init__(self,Args):
                #Input argumenets
		self.Args=Args

		if self.Args.OutputFileName is None:
			self.Args.OutputFileName=self.Args.InputFileName.replace(".vtu","_Plots.dat")

		self.MBF_data={"LAD":[], "Diag1":[], "Diag2":[], "Intermedius":[], "LCx":[], "OM":[], "L_PL":[], "R_PL":[], "PDA":[]}
	
		#The LV filename that contains MBF
		self.VentricleFileName=self.Args.InputFileName

                #Coronary Centerline Files
		self.CenterlineFileNamesPaths=sorted(glob("%s/*.vtp"%self.Args.CenterlinesFolder))
		self.CenterlineFileNames=[filename.split("/")[-1] for filename in self.CenterlineFileNamesPaths]

                ##Separate the Coronary Territories into Left and Right side - Only when using mesh-complete folder
		#self.CenterlineFileNamesPaths=[]
		#self.CenterlineFileNames=[]
		#for i in range(len(CenterlineFileNames)):
		#	CLName_=CenterlineFileNames[i]
		#	if CLName_.find("RCA")>=0 or CLName_.find("LCx")>=0 or CLName_.find("LAD")>=0 or CLName_.find("Ramus")>=0 or CLName_.find("Diag")>=0:
		#		self.CenterlineFileNames.append(CenterlineFileNames[i])
		#		self.CenterlineFileNamesPaths.append(CenterlineFileNamesPaths[i])
		#		print ("------ Loading Centerline: %s"%CenterlineFileNamesPaths[i])"""

	def Main(self):
		#Read the LV mesh file that contains the territory mappings and MBF
		print ("--- Reading Left Ventricle MBF Territories: %s"%self.VentricleFileName)
                
		if   self.Args.InputFileName[-4:]==".vtk":
			VentricleMesh=ReadVTKFile(self.Args.InputFileName)   #Read a VTK volume stack
			VentricleMesh=ThresholdByUpper(VentricleMesh,self.Args.ArrayName,1) #Convert to an unstructured grid 
		elif self.Args.InputFileName[-4:]==".vtu":
			VentricleMesh=ReadVTUFile(self.Args.InputFileName) #Read a VTU unstructured grid
		else:
			print ("The extension %s is not valid for volume"%self.Args.InputFileName[-4:])
			print ("Exiting...")
			exit(1)

		#Convert VTK data into Numpy
		DataArray=vtk_to_numpy(VentricleMesh.GetPointData().GetArray(self.Args.ArrayName))
		
		#Compute Data Statistics
		print ("--- Computing Data Statistics")
		self.Mean=np.average(DataArray)
		self.Stdev=np.std(DataArray)
		self.Mode=float(stats.mode(DataArray)[0]) 
		self.Perc50=float(np.percentile(DataArray,50))
		self.Perc75=float(np.percentile(DataArray,75))
		self.Koen80=float(0.8*self.Perc75)  

		print ("------ Mean                : %.05f"%self.Mean)
		print ("------ Stdev               : %.05f"%self.Stdev)
		print ("------ Mode                : %.05f"%self.Mode)
		print ("------ Median              : %.05f"%self.Perc50)
		print ("------ 75th Perc           : %.05f"%self.Perc75)
		print ("------ 80%% of 75th Perc   : %.05f"%self.Koen80)
		
		#Compute the vessel-specific MBF
		print ("--- Computing vessel-specific Territory Statistics")
		MBF_Raw,MBF_Mean,MBF_Mode,MBF_Perc50,MBF_Perc75,MBF_Koen80=self.compute_mbf(VentricleMesh)
	
	
		#Write the tecplot file 
		#self.write_MBF_data(MBF_Array1)
		self.write_tecplot(MBF_Raw,MBF_Mean,MBF_Mode,MBF_Perc50,MBF_Perc75,MBF_Koen80)

	def compute_mbf(self,VentricleMesh):
                #Separate the vessels into various territories
		MBF_Raw_Array=copy.deepcopy(self.MBF_data)
		MBF_Mean_Array=copy.deepcopy(self.MBF_data)
		MBF_Mode_Array=copy.deepcopy(self.MBF_data)
		MBF_Perc50_Array=copy.deepcopy(self.MBF_data)
		MBF_Perc75_Array=copy.deepcopy(self.MBF_data)
		MBF_Koen80_Array=copy.deepcopy(self.MBF_data)


		#Get the number of points in the array
		N=VentricleMesh.GetNumberOfPoints()

                #Loop over all points and append to the correct array
		for i in range(N):
			Value_=VentricleMesh.GetPointData().GetArray(self.Args.ArrayName).GetValue(i)
			territory_idx_=VentricleMesh.GetPointData().GetArray("TerritoryMaps").GetValue(i)
			filename_=self.CenterlineFileNames[territory_idx_]
			
			for Key_ in self.MBF_data:
				if filename_.find(Key_)>=0:
					MBF_Raw_Array[Key_].append(Value_)
					MBF_Mean_Array[Key_].append(Value_/self.Mean)
					MBF_Mode_Array[Key_].append(Value_/self.Mode)
					MBF_Perc50_Array[Key_].append(Value_/self.Perc50)
					MBF_Perc75_Array[Key_].append(Value_/self.Perc75)
					MBF_Koen80_Array[Key_].append(Value_/self.Koen80)
		
		return MBF_Raw_Array, MBF_Mean_Array,MBF_Mode_Array,MBF_Perc50_Array,MBF_Perc75_Array,MBF_Koen80_Array

       
	def write_tecplot(self,MBF_Raw,MBF_Mean,MBF_Mode,MBF_Perc50,MBF_Perc75,MBF_Koen80):
		outfile=open(self.Args.OutputFileName,'w')
		outfile.write('TITLE="MBF(mL/min/100g)"\n')
		outfile.write('VARIABLES = "Region","Raw", "NormMean", "NormMode", "Norm50Perc", "Norm75Perc", "NormKoen80"\n')

		count=self.Args.XaxisOffset

                #Make Whisker plots for each key
		for key in self.MBF_data:
			print ("------ Looping over Artery: %s"%key)
			if sum(MBF_Mean[key])==0.0: 
				print ("--------- Above Artery Not Found. Skipping")
				count+=1
				continue
			
			Mean_NormMean  =np.mean(MBF_Mean[key])
			Mean_NormMode  =np.mean(MBF_Mode[key])
			Mean_Norm50Perc=np.mean(MBF_Perc50[key])
			Mean_Norm75Perc=np.mean(MBF_Perc75[key])
			Mean_NormKoen80=np.mean(MBF_Koen80[key])

			Max_Raw       =np.max(MBF_Raw[key])
			Max_NormMean  =np.max(MBF_Mean[key])
			Max_NormMode  =np.max(MBF_Mode[key])
			Max_Norm50Perc=np.max(MBF_Perc50[key])
			Max_Norm75Perc=np.max(MBF_Perc75[key])
			Max_NormKoen80=np.max(MBF_Koen80[key])

			Lower_Raw       =np.percentile(MBF_Raw[key],25)
			Lower_NormMean  =np.percentile(MBF_Mean[key],25)
			Lower_NormMode  =np.percentile(MBF_Mode[key],25)
			Lower_Norm50Perc=np.percentile(MBF_Perc50[key],25)
			Lower_Norm75Perc=np.percentile(MBF_Perc75[key],25)
			Lower_NormKoen80=np.percentile(MBF_Koen80[key],25)

                        
			Mid_Raw       =np.percentile(MBF_Raw[key],50)
			Mid_NormMean  =np.percentile(MBF_Mean[key],50)
			Mid_NormMode  =np.percentile(MBF_Mode[key],50)
			Mid_Norm50Perc=np.percentile(MBF_Perc50[key],50)
			Mid_Norm75Perc=np.percentile(MBF_Perc75[key],50)
			Mid_NormKoen80=np.percentile(MBF_Koen80[key],50)

			Upper_Raw       =np.percentile(MBF_Raw[key],75)
			Upper_NormMean  =np.percentile(MBF_Mean[key],75)
			Upper_NormMode  =np.percentile(MBF_Mode[key],75)
			Upper_Norm50Perc=np.percentile(MBF_Perc50[key],75)
			Upper_Norm75Perc=np.percentile(MBF_Perc75[key],75)
			Upper_NormKoen80=np.percentile(MBF_Koen80[key],75)

			Min_Raw       =np.min(MBF_Raw[key])
			Min_NormMean  =np.min(MBF_Mean[key])
			Min_NormMode  =np.min(MBF_Mode[key])
			Min_Norm50Perc=np.min(MBF_Perc50[key])
			Min_Norm75Perc=np.min(MBF_Perc75[key])
			Min_NormKoen80=np.min(MBF_Koen80[key])
			
			outfile.write('Zone T= "%s_%s", I=5, F=POINT\n'%(key,"Markers"))
			outfile.write("%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n"%(count,Min_Raw,Min_NormMean,Min_NormMode,Min_Norm50Perc,Min_Norm75Perc,Min_NormKoen80))
			outfile.write("%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n"%(count,Lower_Raw,Lower_NormMean,Lower_NormMode,Lower_Norm50Perc,Lower_Norm75Perc,Lower_NormKoen80))
			outfile.write("%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n"%(count,Mid_Raw,Mid_NormMean,Mid_NormMode,Mid_Norm50Perc,Mid_Norm75Perc,Mid_NormKoen80))
			outfile.write("%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n"%(count,Upper_Raw,Upper_NormMean,Upper_NormMode,Upper_Norm50Perc,Upper_Norm75Perc,Upper_NormKoen80))
			outfile.write("%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n"%(count,Max_Raw,Max_NormMean,Max_NormMode,Max_Norm50Perc,Max_Norm75Perc,Max_NormKoen80))
			

			#Make a box around the 25-75 percentile for Raw MBF
			outfile.write('Zone T= "%s_%s", I=5, F=POINT\n'%(key,"Box"))
			outfile.write('%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n'%(count-0.25,Lower_Raw,Lower_NormMean,Lower_NormMode,Lower_Norm50Perc,Lower_Norm75Perc,Lower_NormKoen80))
			outfile.write('%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n'%(count+0.25,Lower_Raw,Lower_NormMean,Lower_NormMode,Lower_Norm50Perc,Lower_Norm75Perc,Lower_NormKoen80))
			outfile.write('%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n'%(count+0.25,Upper_Raw,Upper_NormMean,Upper_NormMode,Upper_Norm50Perc,Upper_Norm75Perc,Upper_NormKoen80))
			outfile.write('%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n'%(count-0.25,Upper_Raw,Upper_NormMean,Upper_NormMode,Upper_Norm50Perc,Upper_Norm75Perc,Upper_NormKoen80))
			outfile.write('%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n'%(count-0.25,Lower_Raw,Lower_NormMean,Lower_NormMode,Lower_Norm50Perc,Lower_Norm75Perc,Lower_NormKoen80))
			
			outfile.write('Zone T= "%s_%s", I=2, F=POINT\n'%(key,"Median_Line"))
			outfile.write('%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n'%(count-0.25,Mid_Raw,Mid_NormMean,Mid_NormMode,Mid_Norm50Perc,Mid_Norm75Perc,Mid_NormKoen80))
			outfile.write('%.02f %.05f %.05f %.05f %.05f %.05f %.05f\n\n'%(count+0.25,Mid_Raw,Mid_NormMean,Mid_NormMode,Mid_Norm50Perc,Mid_Norm75Perc,Mid_NormKoen80))
                       

			count+=1
			

		outfile.close()
	

	# This function is to compare pre- and post-cabg data 
	def write_MBF_data(self,MBF_Array1,MBF_Array2,MBF_Array3,MBF_Array4):
                
		outfile=open(self.Args.OutputFileName,'w')
		outfile.write('TITLE="MBF(mL/min/100g)"\n')
		outfile.write('VARIABLES = "Region","%s", "%s", "%s", "%s"\n'%(self.Args.Array1,self.Args.Array2,self.Args.Array3,self.Args.Array4))
                
		if self.Args.PostCABG==0: 
			count=1
			tag="pre"
		else: 
			count=1.25
			tag="post"

		for key in self.MBF_data:
			outfile.write('Zone T= "%s_%s", I=3, F=POINT\n'%(key,tag))
			mean_Array1=np.mean(MBF_Array1[key])
			mean_Array2=np.mean(MBF_Array2[key])
			mean_Array3=np.mean(MBF_Array3[key])
			mean_Array4=np.mean(MBF_Array4[key])

			stdev_Array1=np.std(MBF_Array1[key])
			stdev_Array2=np.std(MBF_Array2[key])
			stdev_Array3=np.std(MBF_Array3[key])
			stdev_Array4=np.std(MBF_Array4[key])

			outfile.write("%.02f %.05f %.05f %.05f %.05f\n"%(count,mean_Array1+stdev_Array1,mean_Array2+stdev_Array2,mean_Array3+stdev_Array3,mean_Array4+stdev_Array4))
			outfile.write("%.02f %.05f %.05f %.05f %.05f\n"%(count,mean_Array1,mean_Array2,mean_Array3,mean_Array4))
			outfile.write("%.02f %.05f %.05f %.05f %.05f\n"%(count,mean_Array1-stdev_Array1,mean_Array2-stdev_Array2,mean_Array3-stdev_Array3,mean_Array4-stdev_Array4))
			count+=1
		outfile.close()

	def Read_Vtu(self,filename):
		Volume=vtk.vtkXMLUnstructuredGridReader()
		Volume.SetFileName(filename)
		Volume.Update()
		Volume=Volume.GetOutput()
		N_p_vol=Volume.GetNumberOfPoints()

		#This options strips a single node from the surface
		#to avoid artifacts close to the surface of endo- or
		#epicardium
		if self.Args.StripEdges!=0:
			#Get the surface from the volume
			Surface=vtk.vtkDataSetSurfaceFilter()
			Surface.SetInputData(Volume)
			Surface.Update()
			Surface=Surface.GetOutput()
			
			#Get all of the surface nodes:
			N_p_surf=Surface.GetNumberOfPoints()
			Surface_Coords_x=np.zeros(N_p_surf)
			Surface_Coords_y=np.zeros(N_p_surf)
			Surface_Coords_z=np.zeros(N_p_surf)
			for i in range(N_p_surf):
				surface_coord_=Surface.GetPoint(i)[:]
				Surface_Coords_x[i]=surface_coord_[0]
				Surface_Coords_y[i]=surface_coord_[1]
				Surface_Coords_z[i]=surface_coord_[2]
			
			Surface_Coords_list_x=np.ndarray.tolist(Surface_Coords_x)
			Surface_Coords_list_y=np.ndarray.tolist(Surface_Coords_y)
			Surface_Coords_list_z=np.ndarray.tolist(Surface_Coords_z)

			#Tag all of the surface coordinates on the volume
			#Make the values on the Surface equal to 0
			Volume_Surface_tags=[]
			for i in range(0,N_p_vol):
				volume_coord_=Volume.GetPoint(i)[:]
				if (volume_coord_[0] in Surface_Coords_list_x) and (volume_coord_[1] in Surface_Coords_list_y) and (volume_coord_[2] in Surface_Coords_list_z):
					Volume_Surface_tags.append(True)
					Surface_Coords_list_x.remove(volume_coord_[0])
					Surface_Coords_list_y.remove(volume_coord_[1])
					Surface_Coords_list_z.remove(volume_coord_[2])
					Volume_Surface_tags.append(True)
				else:
					Volume_Surface_tags.append(False)
			

			exit(1)
		return reader

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will generate tecplot files from perfusion territory surface files.")

        #Input filename of the perfusion map
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=False, dest="InputFileName", default="Results/MBF_ProjectedToSurface_Territories.vtp",help="Volumetric Mesh that contains the Myocardial Blood Flow Data and Territory Maps")

        #Input folder that contains the coronary centerlines
	parser.add_argument('-CenterlinesFolder', '--CenterlinesFolder', type=str, required=False, default="Results/Centerlines", dest="CenterlinesFolder",help="Folder that contains the centerline files")

	#Option to strip off one endo- and one epi-cardial element
	#parser.add_argument('-stripedges', '--StripEdges', type=int, required=False, default=0, dest="StripEdges", help="The option will strip MBF values from exterior and interior of the myocardium")
	#parser.add_argument('-edgespacing', '--EdgeSpacing', type=float, required=False, default=0, dest="EdgeSpacing", help="The distance from interior wall to remove MBF values")	

	##Pre- or Post-CABG
	#parser.add_argument('-postCABG', '-PostCABG', type=int, required=True, default=0, dest="PostCABG", help="0 for pre-cabg (default) and 1 for post-cabg")

        #Output argumenets
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output filename to store the perfusion plots and histograms")
	
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, default="scalars", dest="ArrayName",help="The array name contain the MBF values. Default is 'scalars'")
	
	parser.add_argument('-XaxisOffset', '--XaxisOffset', type=float, required=False, default=1.0, dest="XaxisOffset",help="Offset for the x axis when ploting in tecplot. helpful for Post-Cabag")
	
	args=parser.parse_args()
	
	#if args.StripEdges!=0 and args.EdgeSpacing==0:
	#	parser.error('EdgeSpacing not provided while StripeEdge is turned on')

	CoronaryPerfusionPlot(args).Main()	
