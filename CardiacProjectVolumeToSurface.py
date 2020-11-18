#This script will project quantify the volumetric quantities
#and project this onto the 2D surface
import os
import sys
from glob import glob
import vtk
import numpy as np
import argparse
from vtk.util.numpy_support import vtk_to_numpy as npvtk
from vtk.util.numpy_support import numpy_to_vtk
from time import time as TIME
from PrincipleComponentAnalysis import PRINCIPLE_COMPONENT_ANALYSIS

class CardiacProjectVolumeToSurface():
	def __init__(self,Args):
		#Infolder with all of the data
		self.Args=Args

	def main(self):
		#The input arguments
		Args=self.Args
		
		print ("Reading: ",Args.InputFileName)
		Volume=self.READ_VTU(Args.InputFileName)
		
		print ("--- Getting the Surface")
		Surface=self.GET_SURFACE(Volume)
	
		if Args.Smoothing==1: 	
			print ("--- Smoothing the Surface")
			Surface=self.SURFACE_SMOOTHING(Surface,Args.Iterations,Args.PassBand,method="Laplace")
		print ("--- Computing Princple Component Analysis and writing file")
		Centroid,Norm1,Norm2,Coord_apex,Size=PRINCIPLE_COMPONENT_ANALYSIS(Volume).main()

		print ("--- Getting Coordinates and Surface Normals")
		Coords,SurfaceNormals=self.GET_COORDINATES(Surface,Centroid)

		print ("--- Create the Slices of the Myocardium")
		Slices,Centroid_slices=self.MYOCARDIUM_SLICES(Volume,Centroid,Norm1,Size)

                #Tag the outer section of the ventricle (0=Outer,  1==Inner)
		print ("--- Extracting outer surface of the ventricle (i.e. epicardium)")
		Surface=self.TAG_OUTER_SURFACE(Surface,Coords,SurfaceNormals,Size)

		#Draw a Line from each Node to Centerline and Get Average Quantities
		print ("--- Computing the Average MBF on the outer wall.")
		Line_avg=self.GET_LINE_AVERAGE(Surface,Slices,Centroid,Coords,Size,Norm1,Norm2,Centroid_slices)
		#Store the line average data back to surface
		Surface_avg=self.SET_LINE_DATA(Surface,Line_avg)

		#Store the Data to the output
		self.WRITE_VTP(Surface_avg,Args.OutputFileName)	
	
	#Set the Line average data to Surface
	def SET_LINE_DATA(self,Surface,Line_avg):
		Surface_avg=numpy_to_vtk(Line_avg,deep=True)
		Surface_avg.SetName("%s_Projected"%self.Args.ArrayName)
		Surface.GetPointData().AddArray(Surface_avg)
		Surface.Modified()
		return Surface

	#Draw a Line from each Node to Centerline
	def GET_LINE_AVERAGE(self,Surface,Slices,Centroid,Coords,Size,Norm1,Norm2,Centroid_slices):
		Point_base=Centroid+2*Size*Norm1
		Point_apex=Centroid-2*Size*Norm1
		Npts=len(Coords)
		
		#Data array to store the average value along the line
		Average_Data=np.zeros(Npts)	
		progress_old=-1
		Time_array=[]
		for i in range(Npts):
			progress_=self.PRINT_PROGRESS(i,Npts,progress_old)
			progress_old=progress_
			if Surface.GetPointData().GetArray("Tags").GetValue(i)==1:
				continue
			else:
				progress_=self.PRINT_PROGRESS(i,Npts,progress_old)
				progress_old=progress_
				#Distance from Pn to Apex
				D_hyp=np.linalg.norm(Coords[i]-Point_apex)
				#Distance from Pn to Centerline (vtk function returns distance squared)
				D_adj=(vtk.vtkLine().DistanceToLine(Coords[i], Point_apex, Point_base))**0.5
				#Projected Distance along the Centerline
				D_opp = np.power(np.power(D_hyp,2)-np.power(D_adj,2),0.5)
			
				#Now compute the projected point on the centerline
				pSource_=Coords[i] #Point on the outer wall
				pTarget_=Point_apex+D_opp*Norm1 #Point on the CL
			
				#Clip the surface to on the plane of the outer point
				#slice_=self.CLIPPED_SLICE(pTarget_,Norm1,Volume) #Slice the plane
				n_=int(D_opp/self.Args.Resolution)
				slice_=Slices[n_]
				centroid_slice_=Centroid_slices[n_]
				line_,time_=self.CLIPPED_LINE(slice_,pSource_,pTarget_,Norm1,centroid_slice_) 
				Time_array.append(time_)
				#print (np.average(Time_array))
				#Compute the average along the line
				n_points_line_=line_.GetNumberOfPoints()
				if n_points_line_>0:
					line_av_=np.zeros(n_points_line_)
					for j in range(n_points_line_):
						line_av_[j]=line_.GetPointData().GetArray(self.Args.ArrayName).GetValue(j)
					Average_Data[i]=np.mean(line_av_)
		return Average_Data

        
	def CLIPPED_LINE(self,Slice,pSource,pTarget,Norm1,centroid_slice):
		time_0=TIME()
                #Get the two in-plane normals
		Norm2_slice=(pSource-centroid_slice)/np.linalg.norm(pSource-centroid_slice)
		Norm3_slice=np.cross(Norm1,Norm2_slice)
		#Generate the two planes
		plane_N2=vtk.vtkPlane()
		plane_N2.SetOrigin(centroid_slice)
		plane_N2.SetNormal(Norm2_slice)
		plane_N3=vtk.vtkPlane()
		plane_N3.SetOrigin(centroid_slice)
		plane_N3.SetNormal(Norm3_slice)
		#Clip the plane to get a line across the diameter
		Line =vtk.vtkCutter()
		Line.GenerateTrianglesOff()
		Line.SetCutFunction(plane_N3)
		Line.SetInputData(Slice)
		Line.Update()
		#Separate the line into only one quarter (i.e. half the line)
		Line1=vtk.vtkClipPolyData()
		Line1.SetClipFunction(plane_N2)
		Line1.SetInputData(Line.GetOutput())
		Line1.Update()
		Line1_data=Line1.GetOutput()

		del plane_N2,plane_N3,Line,Line1
		time_=TIME()-time_0
		return Line1_data,time_

        
	def MYOCARDIUM_SLICES(self,Volume,Centroid,Norm1,Size):
		#Get the base and apex coordinate
		Coord_base=Centroid+2*Size*Norm1
		Coord_apex=Centroid-2*Size*Norm1
		#Create slice along the length of the appex
		Npts=int((4*Size)/self.Args.Resolution)
		#Loop over the the myocardium and compute slices        
		Slices=[]
		Centroid_slices=[]
		progress_old=-1
		for i in range(Npts):
			progress_=self.PRINT_PROGRESS(i,Npts,progress_old)
			progress_old=progress_
			coord_=Coord_apex+Norm1*i*self.Args.Resolution
			slice_=self.CLIPPED_SLICE(coord_,Norm1,Volume)
			if slice_.GetNumberOfPoints()==0:
				Centroid_slices.append(coord_)
			else: Centroid_slices.append(self.GET_CENTROID(slice_))
			Slices.append(slice_)
		return Slices,Centroid_slices

	#Get the outer surface of the mesh
	def TAG_OUTER_SURFACE(self,Surface,Coords,SurfaceNormals,Size):
		Npts=len(Coords)
		#Get the P1 pointing away from the Surface Nodes
		Coords_p1=Coords+SurfaceNormals*2*Size
		#Create an OBB tree and cast Rays	
		obbTree = vtk.vtkOBBTree()
		obbTree.SetDataSet(Surface)
		obbTree.BuildLocator()
		pointsVTKintersection = vtk.vtkPoints()
		Surface_tags=np.zeros(Npts)
		for i in range(Npts):
			pSource=Coords[i]
			pTarget=Coords_p1[i]
			code = obbTree.IntersectWithLine(pSource, pTarget, pointsVTKintersection, None)
			X=pointsVTKintersection.GetData().GetNumberOfTuples()
			if X>1: 
				Surface_tags[i]=1	
				SurfaceNormals[i]=SurfaceNormals[i]*-1	
		
		#Store the data in Surface array
		#Normals
		Normals_vtk=numpy_to_vtk(SurfaceNormals, deep=True)
		Normals_vtk.SetName("Normals")
		Surface.GetPointData().AddArray(Normals_vtk)
		#Tags for out or inner surface
		Surface_tags_vtk=numpy_to_vtk(Surface_tags,deep=True)
		Surface_tags_vtk.SetName("Tags")
		Surface.GetPointData().AddArray(Surface_tags_vtk)	
		Surface.Modified()

		return Surface

	#Smooth the surface
	def SURFACE_SMOOTHING(self,Surface,Nits,PB_value,method="Taubin"):
		if method=="Taubin":
			smoothingFilter = vtk.vtkWindowedSincPolyDataFilter()
			smoothingFilter.SetInputData(Surface)
			smoothingFilter.SetNumberOfIterations(Nits)
			smoothingFilter.SetPassBand(PB_value)
			smoothingFilter.SetBoundarySmoothing(True)
			smoothingFilter.Update()
			return smoothingFilter.GetOutput() 
		elif method=="Laplace":
			smoothingFilter = vtk.vtkSmoothPolyDataFilter()
			smoothingFilter.SetInputData(Surface)
			smoothingFilter.SetNumberOfIterations(Nits)
			smoothingFilter.SetRelaxationFactor(PB_value)
			smoothingFilter.Update()
			return smoothingFilter.GetOutput()
		else:
			print ("Error. The smoothing filter was not found")
			exit(1)

	#Get the Surface of the Mesh
	def GET_SURFACE(self,Volume):
		Surface=vtk.vtkDataSetSurfaceFilter()
		Surface.SetInputData(Volume)
		Surface.Update()
		return Surface.GetOutput()

        #Read the VTU file
	def READ_VTU(self,filename):
		Volume=vtk.vtkXMLUnstructuredGridReader()
		Volume.SetFileName(filename)
		Volume.Update()
		return Volume.GetOutput()

	def WRITE_VTP(self,Surface,Filename):
		writer=vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(Filename)
		writer.SetInputData(Surface)
		writer.Update()
        
        #Get Centroid of the VTK dataset
	def GET_CENTROID(self,Surface):
		Centroid=vtk.vtkCenterOfMass()
		Centroid.SetInputData(Surface)
		Centroid.SetUseScalarsAsWeights(False)
		Centroid.Update()
		return Centroid.GetCenter()

	def GET_COORDINATES(self,Surface,Centroid):
                #Loop over all the points
		Npts=Surface.GetNumberOfPoints()
		Coords=np.zeros(shape=(Npts,3))
		for i in range(Npts): Coords[i,:]=Surface.GetPoints().GetPoint(i)[:]
		#Get the vector from Centroid to Surface Node and normalize them
		Normals=Coords-Centroid
		for i in range(Npts): Normals[i]=Normals[i]/np.linalg.norm(Normals[i])
		return Coords,Normals

        
	def CLIPPED_SLICE(self,Origin,Norm,Volume):
		plane=vtk.vtkPlane()
		plane.SetOrigin(Origin)
		plane.SetNormal(Norm)
		Slice=vtk.vtkCutter()
		Slice.GenerateTrianglesOff()
		Slice.SetCutFunction(plane)
		Slice.SetInputData(Volume)
		Slice.Update()
		del plane
		return Slice.GetOutput()
	
	def GET_THRESHOLD_SURFACE(self,Surface,array_name,value):
		thresh = vtk.vtkThreshold()
		thresh.SetInputData(Surface)
		thresh.ThresholdByUpper(value)
		thresh.SetInputArrayToProcess(0, 0, 0,"vtkDataObject::FIELD_ASSOCIATION_POINTS",array_name)
		thresh.Update()
		vtk_grid = thresh.GetOutput()	
		return vtk_grid
	
	#Print the progress of the loop
	def PRINT_PROGRESS(self,i,N,progress_old):
		progress_=(int((float(i)/N*100+0.5)))
		if progress_%10==0 and progress_%10!=progress_old: print ("    Progress: %d%%"%progress_)
		#progress_=int((float(i)/N)*1000)
		#if progress_%100==0: print ("    Progress: %d%%"%int(progress_/10.)),
		return progress_%10	
if __name__=='__main__':
	#Description
	parser = argparse.ArgumentParser(description="This script will project the volumetric Myocardial Blood Flow data onto the surface of the ventricle.")
	
	#Input filename of the perfusion map
	parser.add_argument('-ifile', '--InputFileName', type=str, required=True, dest="InputFileName",help="Volumetric Mesh that contains the Myocardial Blood Flow Data ")
	
	#Resolution of the Averaging Procedure
	parser.add_argument('-resolution', '--Resolution', type=float, required=False,default=0.02, dest="Resolution",help="The resolution of the averaging procedure (default=0.02)")
	
	#Surface Smoothing
	parser.add_argument('-smoothing', '--Smoothing', type=int, required=False,default=1, dest="Smoothing",help="Smooth the surface (default=True)")
	parser.add_argument('-passband', '--PassBand', type=float, required=False,default=0.2, dest="PassBand",help="Passband value for the smoothing (default=0.2)")
	parser.add_argument('-iterations', '--Iterations', type=int, required=False,default=500, dest="Iterations",help="Iterations for the smoothing (default=500)")
	
	#Array Name of the Data
	parser.add_argument('-arrayname', '--ArrayName', type=str, required=False,default="ImageScalars", dest="ArrayName",help="The name of the array containing the MBF values")
	
	#Output argumenets
	parser.add_argument('-ofile', '--OutputFileName', type=str, required=True, dest="OutputFileName",help="The output filename of the projected surface")

        
	args=parser.parse_args()
	CardiacProjectVolumeToSurface(args).main()
