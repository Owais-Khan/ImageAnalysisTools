import numpy
import argparse
import parser
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from utilities import *

class ImageAnalysisVTKToRaw():
	def __init__(self,Args):
		self.Args=Args
	def Main(self):
		#Load the VTI File 
		print ("Reading the File: %s"%self.Args.InputFileName)
		ImageDataVTI=ReadVTIFile(self.Args.InputFileName)
		
		#Get the number of points
		Npts=ImageDataVTI.GetNumberOfPoints()
		print ("The number of data points are: %d"%Npts)
		
		#Get the shape of the image stack
		Dims=ImageDataVTI.GetDimensions()
		print ("The dimensions of the image are:",Dims)

		#Get the origin of the image
		Origin=ImageDataVTI.GetOrigin()
		print ("The origin of the image is: ",Origin)

		#Get the spacing of the image
		Spacing=ImageDataVTI.GetSpacing()
		print ("The spacing of the image is: ",Spacing)
			
	
		Labels=vtk_to_numpy(ImageDataVTI.GetPointData().GetArray("Labels"))
		ImageValues=vtk_to_numpy(ImageDataVTI.GetPointData().GetArray("Scalars_"))
		
		#Reshape the Labels	
		if self.Args.ReshapeArray==1:
			print ("Reshaping the arrays to image dimensions...")
			Labels=np.reshape(Labels, Dims)
			ImageValues=np.reshape(ImageValues,Dims)	

		#Save the file in .raw format
		outfile=self.Args.InputFileName.replace(".vti",".raw")
		print ("Writing Image Scalars to .raw file: %s"%outfile)
		ImageValues.astype('int16').tofile(outfile)
		outfile_labels=outfile.replace(".raw","_Labels.raw")
		print ("Writing Image Scalars to .raw file: %s"%outfile_labels)
		Labels.astype("int16").tofile(outfile_labels)


if __name__=="__main__":
        parser = argparse.ArgumentParser(description="This script will convert the vtk image stack into 2D slices written in .raw format.")

        parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The file name containing the .vtk data that needs to be written to .raw format")
       
        parser.add_argument('-ReshapeArray', '--ReshapeArray', type=int, required=False,default=0, dest="ReshapeArray",help="1=Reshape the array to original image stack shape, 0 (default)=Keep a 0D array of size=IxJxK")

        args=parser.parse_args()
        ImageAnalysisVTKToRaw(args).Main()
 
