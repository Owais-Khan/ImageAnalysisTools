#This script will seperate the images into
#individual cycles for perfusion cycle
#and extract the cycle with highest 
#contrast intensity.
#The cycle will be placed under ./casename/dicom_perfusion_cycle_X"
from __future__ import print_function
import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import numpy as np
import vtk
import argparse
from matplotlib import pyplot as plt
class CardiacSeparatePerfusionCycles():
	def __init__(self,Args):
		#Store the input arguments
		self.Args=Args

		#Seperate the files into individual folders
		filenames=glob("%s/*.dcm"%self.Args.InputFolderName)
		try:
			filenames=sorted(filenames, key=lambda filename: int(filename.split(".")[4]))
		except: 
			filenames=sorted(filenames)
		N_image_cyc=int((len(filenames))/self.Args.NumberOfCycles)
		count=0
		for i in range(0,len(filenames),N_image_cyc):
			directoryname="%s/cycle_%d/"%(self.Args.InputFolderName,count)
			if len(glob(directoryname))==0:
				os.system("mkdir %s/cycle_%d/"%(self.Args.InputFolderName,count))
			#Move the files into separate folder
			for j in range(i,i+N_image_cyc):
				os.system("cp %s %s/cycle_%d/"%(filenames[j],self.Args.InputFolderName,count))
				
			count+=1

		print ("/n")
		print ("The number of cycles are: %d"%self.Args.NumberOfCycles)
		print ("The total number of dicom files are: %d"%len(filenames))
		print ("Number of files per cycle are: %d"%N_image_cyc)
		

	def main(self):				
		#Now get the pixel intensity of each cycle
		#to get the most IV injection curve
		mean_pix_intensity=[]
		for i in range(self.Args.NumberOfCycles):
			print ("--- Looping over cycle: %s"%i)
			#Get Numpy array from the pixel intensity
			PathDicom = "%s/cycle_%s/"%(self.Args.InputFolderName,i)
			reader = vtk.vtkDICOMImageReader()
			reader.SetDirectoryName(PathDicom)
			reader.Update()
			imageData=reader.GetOutput()
			ArrayData=npvtk(imageData.GetPointData().GetArray(0))
			mean_pix_intensity.append(np.mean(ArrayData[ArrayData>0]))
			del reader,imageData,ArrayData

		#Copy the cycle in the main folder that has 
		#highest constrast intensity
		N_max_cyc=mean_pix_intensity.index(max(mean_pix_intensity))
		os.system("mv %s/cycle_%d %s/cycle_%s_MaxIntensity"%(self.Args.InputFolderName,N_max_cyc,self.Args.InputFolderName,N_max_cyc))
		
		print ("Cycle with max intensity is: %s/cycle_%s_MaxIntensity"%(self.Args.InputFolderName,N_max_cyc))

		#Get the Mean Hounsfield unit for Diastolic CTA
		if self.Args.InputFolderName2 is not None:
			PathDicom="%s"%(self.Args.InputFolderName2)
			reader = vtk.vtkDICOMImageReader()
			reader.SetDirectoryName(PathDicom)
			reader.Update()
			imageData=reader.GetOutput()
			ArrayData=npvtk(imageData.GetPointData().GetArray(0))
			CTA_dias_pix_intensity=np.mean(ArrayData[ArrayData>0])
			del reader,imageData,ArrayData
	
			#Output the folder that has the mean Diastolic Pixel Intensity
			N_dias_cyc=(np.abs(mean_pix_intensity - CTA_dias_pix_intensity)).argmin()
			os.system("mv %s/cycle_%d %s/cycle_%s_AnatomicIntensity"%(self.Args.InputFolderName,N_dias_cyc,self.Args.InputFolderName,N_dias_cyc))
			print ("Cycle with intensity same as anatomic scan is: %s/cycle_%s_AnatomicIntensity"%(self.Args.InputFolderName,N_dias_cyc))
		print ("Contrast Intensity of Cycle %d (max) is: %.05f"%(N_max_cyc,mean_pix_intensity[N_max_cyc])) 
		if self.Args.InputFolderName2 is not None:
			print ("Contrast Intensity of Cycle %d (anatomic) is: %.05f"%(N_dias_cyc,CTA_dias_pix_intensity)) 
		#Get the Mean Hounsfield unit for the Systolic CTA
		"""PathDicom="%s/dicom_CTA_systole/"%(self.casename)
		reader = vtk.vtkDICOMImageReader()
		reader.SetDirectoryName(PathDicom)
		reader.Update()
		imageData=reader.GetOutput()
		ArrayData=npvtk(imageData.GetPointData().GetArray(0))
		CTA_sys_pix_intensity=np.mean(ArrayData[ArrayData>0])		
		del reader,imageData,ArrayData"""

		#Plot the data in matplot lib
		print (mean_pix_intensity)
		plt.plot(mean_pix_intensity,'-ko',label="Perfusion")
		if self.Args.InputFolderName2 is not None: plt.axhline(y=CTA_dias_pix_intensity,color='r',label="Diastole")
		#plt.axhline(y=CTA_sys_pix_intensity,color='b',label="Systole")
		plt.legend()
		plt.xlabel("Cycles")
		plt.ylabel("Hounsfield Units")
		plt.show()
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take a dicom folder with N cycles and separate them into individual folder, and also write out the cycle with peak contrast intensity.")

        #Input filename of the perfusion map
	parser.add_argument('-ifolder', '--InputFolderName', type=str, required=True, dest="InputFolderName",help="The name of the folder where all dicom files are stored for N cycles")
	
	#Number of cycles 
	parser.add_argument('-cycles', '--NumberOfCycles', type=int, required=True, dest="NumberOfCycles",help="The number of cycles that are in the dicom folder")
	
	#Diastolic cycle
	parser.add_argument('-ifolder2', '--InputFolderName2', type=str, required=False, dest="InputFolderName2",help="The name of the folder that contains dicom files for anatomic (diastolic) scan")

	#Output folder
	#parser.add_argument('-ofolder', '--OutputFolderName', type=str, required=True, dest="OutputFolderName",help="The name of the folder where the cycle with peak contrast instensity will be stored")
	
	args=parser.parse_args()
	CardiacSeparatePerfusionCycles(args).main()

