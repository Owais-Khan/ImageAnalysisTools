import vtk
#from vmtk import vtkvmtk, vmtkscripts
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
from glob import glob
from scipy.spatial import distance as DISTANCE
from scipy.stats import iqr as IQR
from scipy.stats import kurtosis as KURTOSIS
from scipy.stats import skew as SKEWNESS
from scipy.stats import mode as MODE


############ Read Dicom Folder ############
"""def ReadDicomFiles(FolderName):
	FileList=sorted(glob("%s/*.dcm"%FolderName))
	print ("------ Number of DICOM Files: %d"%len(
	Image=vmtkscripts.vmtkImageReader()
	Image.InputFileName(FileList[0])
	Image.ImageOutputFileName("/Users/mokhan/GoogleDrive/Owais/Research_Postdoc/perfusion_project/Simvascular/CABG1A/Images/abc.vti")
	Image.Update()

#ReadDicomFiles("/Users/mokhan/GoogleDrive/Owais/Research_Postdoc/perfusion_project/Simvascular/CABG1A/Images/CTA")"""

	
	

############ Input/Output ##################
def ReadVTUFile(FileName):
	reader=vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(FileName)
	reader.Update()
	return reader.GetOutput()

def ReadVTKFile(FileName):
	reader = vtk.vtkStructuredPointsReader()
	reader.SetFileName(FileName)
	reader.ReadAllVectorsOn()
	reader.ReadAllScalarsOn()
	reader.Update()
	return reader.GetOutput()

def ReadVTPFile(FileName):
	reader=vtk.vtkXMLPolyDataReader()
	reader.SetFileName(FileName)
	reader.Update()
	return reader.GetOutput()

def ReadVTIFile(FileName):
	reader = vtk.vtkXMLImageDataReader() 
	reader.SetFileName(FileName) 
	reader.Update()
	return reader.GetOutput()

def WriteVTIFile(FileName,Data):
	writer=vtk.vtkXMLImageDataWriter()
	writer.SetFileName(FileName)
	writer.SetInputData(Data)
	writer.Update()

def WriteVTUFile(FileName,Data):
	writer=vtk.vtkXMLUnstructuredGridWriter()
	writer.SetFileName(FileName)
	writer.SetInputData(Data)
	writer.Update()
        
def WriteVTPFile(FileName,Data):
	writer=vtk.vtkXMLPolyDataWriter()
	writer.SetFileName(FileName)
	writer.SetInputData(Data)
	writer.Update()

def WriteNiftiFile(FileName,Data):
	writer=vtk.vtkNIFTIImageWriter()
	writer.SetFileName(FileName)
	writer.SetInputData(Data)
	writer.Update()

############# Mesh Morphing Functions ###############
        #Create a line from apex and centroid of the myocardium
        
def CreateLine(Point1,Point2,Length):
	line0=np.array([Point1[0]-Point2[0],Point1[1]-Point2[1],Point1[2]-Point2[2]])
	line1=-1*line0
	line0=(line0/np.linalg.norm(line0))*(Length/2.)
	line1=(line1/np.linalg.norm(line1))*(Length/2.)
	return line0,line1

def CreatePolyLine(Coords):
	# Create a vtkPoints object and store the points in it
	points = vtk.vtkPoints()
	for i in range(len(Coords)): points.InsertNextPoint(Coords[i])

	#Create a Polyline
	polyLine = vtk.vtkPolyLine()     
	polyLine.GetPointIds().SetNumberOfIds(len(Coords))
	for i in range(len(Coords)): polyLine.GetPointIds().SetId(i, i)

	# Create a cell array to store the lines in and add the lines to it
	cells = vtk.vtkCellArray()
	cells.InsertNextCell(polyLine)

	# Create a polydata to store everything in
	polyData = vtk.vtkPolyData()
    
	# Add the points to the dataset
	polyData.SetPoints(points)

	# Add the lines to the dataset
	polyData.SetLines(cells)
	
	return polyData 

def ClosestPoint(Point, Array):
	dist_2 = np.sum((Array - Point)**2, axis=1)
	return Array[np.argmin(dist_2)],np.argmin(dist_2),min(dist_2)

def FurthestPoint(Point, Array):
        dist_2 = np.sum((Array - Point)**2, axis=1)
        return Array[np.argmax(dist_2)],np.argmax(dist_2)

        
def CutPlane(Volume,Origin,Norm):
	plane=vtk.vtkPlane()
	plane.SetOrigin(Origin)
	plane.SetNormal(Norm)
	Slice=vtk.vtkCutter()
	Slice.GenerateTrianglesOff()
	Slice.SetCutFunction(plane)
	Slice.SetInputData(Volume)
	Slice.Update()
	return Slice.GetOutput()

def CutLine(Slice,Point,Centroid,Norm1):
	#Get the two in-plane normals
	Norm2_slice=(Point-Centroid)/np.linalg.norm(Point-Centroid)
	Norm3_slice=np.cross(Norm1,Norm2_slice)
	
	#Generate the two planes
	plane_N2=vtk.vtkPlane()
	plane_N2.SetOrigin(Centroid)
	plane_N2.SetNormal(Norm2_slice)
	plane_N3=vtk.vtkPlane()
	plane_N3.SetOrigin(Centroid)
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
	return Line1.GetOutput()

def CreateSphere(Coord,Radius):
	Sphere=vtk.vtkSphere()
	Sphere.SetCenter(Coord)
	Sphere.SetRadius(Radius)
	return Sphere

def GetCentroid(Surface):
	Centroid=vtk.vtkCenterOfMass()
	Centroid.SetInputData(Surface)
	Centroid.SetUseScalarsAsWeights(False)
	Centroid.Update()
	return Centroid.GetCenter()


def ExtractSurface(volume):
	#Get the outer surface of the volume
	surface=vtk.vtkDataSetSurfaceFilter()
	surface.SetInputData(volume)
	surface.Update()
	return surface.GetOutput()
        
#Print the progress of the loop
def PrintProgress(i,N,progress_old):
	progress_=(int((float(i)/N*100+0.5)))
	if progress_%10==0 and progress_%10!=progress_old: print ("    Progress: %d%%"%progress_)
	return progress_%10

def TagOuterSurface(Surface):
	#Create an OBB tree and cast Rays       
	obbTree = vtk.vtkOBBTree()
	obbTree.SetDataSet(Surface)
	obbTree.BuildLocator()
	pointsVTKintersection = vtk.vtkPoints()

	#Create an array to store surface tags
	Surface_tags=np.zeros(Surface.GetNumberOfPoints())
       
	#Get Centroid
	Centroid=np.array(GetCentroid(Surface))

	#Loop over all the points. 
	for i in range(Surface.GetNumberOfPoints()):
		pSource=np.array(Surface.GetPoint(i))
		pTarget=pSource+np.array((pSource-Centroid))*5
		code = obbTree.IntersectWithLine(pSource, pTarget, pointsVTKintersection, None)
		X=pointsVTKintersection.GetData().GetNumberOfTuples()
		if X>1: Surface_tags[i]=1
			
	#Store the data in Surface array
	#Tags for out or inner surface
	Surface_tags_vtk=numpy_to_vtk(Surface_tags,deep=True)
	Surface_tags_vtk.SetName("Tags")
	Surface.GetPointData().AddArray(Surface_tags_vtk)
	Surface.Modified()
                
	return Surface

#Smooth Surface
def SurfaceSmoothing(Surface,Nits,PB_value,method="Taubin"):
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

def SurfaceAddArray(Surface,Array,ArrayName):
	SurfaceArray=numpy_to_vtk(Array,deep=True)
	SurfaceArray.SetName(ArrayName)
	Surface.GetPointData().AddArray(SurfaceArray)
	Surface.Modified()
	return Surface

def SurfaceAddCellArray(Surface,Array,ArrayName):
        SurfaceArray=numpy_to_vtk(Array,deep=True)
        SurfaceArray.SetName(ArrayName)
        Surface.GetCellData().AddArray(SurfaceArray)
        Surface.Modified()
        return Surface


def ProjectedPointOnLine(coord_,Centroid,Apex,Norm1):
	#Find the location (coord,distance) on the LV Apex-Base axis
	dist_P_to_line_=np.sqrt(vtk.vtkLine.DistanceToLine(coord_,Centroid,Apex))
	dist_P_to_Apex_=np.power( np.power(coord_[0]-Apex[0],2) + np.power(coord_[1]-Apex[1],2) + np.power(coord_[2]-Apex[2],2),0.5)
	dist_Apex_to_ProjP_=np.power(np.power(dist_P_to_Apex_,2)-np.power(dist_P_to_line_,2),0.5)
	coord_ProjP_=Apex-Norm1*dist_Apex_to_ProjP_

	return coord_ProjP_

def SurfaceNormals(Surface):
	normals = vtk.vtkPolyDataNormals()
	normals.SetInputData(Surface)
	normals.SetFeatureAngle(90)
	normals.AutoOrientNormalsOn()
	#  normals.GetOutput().ReleaseDataFlagOn()
	normals.UpdateInformation()
	normals.Update()
	Surface = normals.GetOutput()

	return Surface

def ThresholdByUpper(Volume,arrayname,value):
	Threshold=vtk.vtkThreshold()
	Threshold.SetInputData(Volume)
	Threshold.ThresholdByUpper(value)
	Threshold.SetInputArrayToProcess(0,0,0,"vtkDataObject::FIELD_ASSOCIATION_POINTS",arrayname)
	Threshold.Update()
	return Threshold.GetOutput()

def ThresholdInBetween(Volume,arrayname,value1,value2):
        Threshold=vtk.vtkThreshold()
        Threshold.SetInputData(Volume)
        Threshold.ThresholdBetween(value1,value2)
        Threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,arrayname)
        Threshold.Update()
        return Threshold.GetOutput()

def ConvertPointsToLine(PointsArray):
        # Create a vtkPoints object and store the points in it
        Points = vtk.vtkPoints()
        for Point_ in PointsArray:
                Points.InsertNextPoint(Point_)

        #Create a Polyline
        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(len(PointsArray))

        for i in range(0, len(PointsArray)):
                polyLine.GetPointIds().SetId(i, i)


        # Create a cell array to store the lines in and add the lines to it
        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polyLine)

        # Create a polydata to store everything in
        polyData = vtk.vtkPolyData()

        # Add the points to the dataset
        polyData.SetPoints(Points)

        # Add the lines to the dataset
        polyData.SetLines(cells)

        return polyData


def Statistics(Volume,ArrayName,NormalizationValue=None):
	statistics={"Volume":None, "Mean":None, "Stdev": None, "MeanNormalized":None, "StdevNormalized":None, "Median": None, "IQR":None, "Mode":None, "75thPerct": None, "Kurtosis": None, "Skewness":None, "Volume":None}

	#Convert VTK array to numpy
	Data_=vtk_to_numpy(Volume.GetPointData().GetArray(ArrayName))
	statistics["Mean"]      =np.mean(Data_)
	statistics["Stdev"]     =np.std(Data_)
	statistics["75thPerct"] =np.percentile(Data_,75)
	statistics["Median"]    =np.median(Data_)
	statistics["IQR"]       =IQR(Data_)
	if NormalizationValue is None:
		statistics["MeanNormalized"]=statistics["Mean"]/statistics["75thPerct"]
		statistics["StdevNormalized"]=statistics["Stdev"]/statistics["75thPerct"]
	else:
		statistics["MeanNormalized"]=statistics["Mean"]/NormalizationValue
		statistics["StdevNormalized"]=statistics["Stdev"]/NormalizationValue

	statistics["Skewness"]   =SKEWNESS(Data_)
	statistics["Kurtosis"]   =KURTOSIS(Data_)
	statistics["Mode"]       =MODE(Data_)[0]	
	

	Mass = vtk.vtkIntegrateAttributes()
	Mass.SetInputData(Volume)
	Mass.Update() 
	MassData=Mass.GetOutput()
	statistics["Volume"]=MassData.GetCellData().GetArray("Volume").GetValue(0)
	return statistics
		

def LargestConnectedRegion(Volume):
	ConnectedVolume=vtk.vtkConnectivityFilter()
	ConnectedVolume.SetInputData(Volume)
	ConnectedVolume.SetExtractionModeToLargestRegion()
	ConnectedVolume.Update()
	ConnectedVolumeData=ConnectedVolume.GetOutput()
	return ConnectedVolumeData
