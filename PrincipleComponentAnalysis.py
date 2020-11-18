import vtk
from sklearn.decomposition import PCA
import numpy as np

class PRINCIPLE_COMPONENT_ANALYSIS():
	def __init__(self,volume):
		self.volume=volume
		#Get the outer surface of the volume
		self.surface=vtk.vtkDataSetSurfaceFilter()
		self.surface.SetInputData(volume)
		self.surface.Update()
		self.volume=self.surface
		
	def main(self):
		volume_data =self.volume.GetOutput()
		surface_data=self.surface.GetOutput()

		#Get the point cloud data
		Centroid,Coords=self.get_point_data(volume_data)
		
		#Perform PCA analysis
		Norm1,Norm2=self.major_axis(Coords)
		Norm1=Norm1*-1
	
		#Get the Apex Point of the left ventricle
		Apex=self.find_apex(Coords,Centroid,Norm1,surface_data)	
	
		#Get the size of the Myocardium (Apex to valve)
		#Size1=np.max(np.power(np.power(Coords[:,0]-Apex[0],2)+np.power(Coords[:,1]-Apex[1],2)+np.power(Coords[:,2]-Apex[2],2),0.5))
		Size=self.get_size(Coords,Apex,Norm1,Norm2)	
		
		#print ("\n")
		#print ("--- LV Morpholocial Analysis")
		print ("----------The Centroid of LV:    ",Centroid)
		print ("----------The Major Axis Norm:   ",Norm1)
		print ("----------The Minor Axis Norm:   ",Norm2)
		print ("----------The Location of Apex:  ",Apex)
		print ("----------The Size of Myocardium:",Size)
		del self.surface,self.volume,volume_data,surface_data,
		return Centroid,Norm1,Norm2,Apex,Size

	def get_size(self,Coords,Apex,Norm1,Norm2):
		#Get the projection of the max distance onto the myocardial centerline	
		distance=np.power(np.power(Coords[:,0]-Apex[0],2)+np.power(Coords[:,1]-Apex[1],2)+np.power(Coords[:,2]-Apex[2],2),0.5)
		max_id=np.argmax(distance)
		coord_max=Coords[max_id,:]
		vector_max=[coord_max[0]-Apex[0],coord_max[1]-Apex[1],coord_max[2]-Apex[2]]
		vector_max_mag=np.linalg.norm(vector_max)
		vector_max_norm=vector_max/vector_max_mag	
		
		#Now get the theta
		cos_theta=np.dot(vector_max_norm,Norm1)
		Size=vector_max_mag*cos_theta
		return Size	
			

	def find_apex(self,Coords,Cent,Norm1,surface_data):
		#Get the min and max and scale the norm
		xmin=min(Coords[:,0]);xmax=max(Coords[:,0])
		ymin=min(Coords[:,1]);ymax=max(Coords[:,1])
		zmin=min(Coords[:,2]);zmax=max(Coords[:,2])
		Length=max(xmax-xmin,ymax-ymin,zmax-zmin)
		Norm1s=Norm1*Length	
		
		#Get the source point
		pSource=[Cent[0]-Norm1s[0],Cent[1]-Norm1s[1],Cent[2]-Norm1s[2]]
		pTarget=[Cent[0]+Norm1s[0],Cent[1]+Norm1s[1],Cent[2]+Norm1s[2]]

		#Create a bounding box of the surface mesh
		obbTree = vtk.vtkOBBTree()
		obbTree.SetDataSet(surface_data)
		obbTree.BuildLocator()
		
		#Create a point array to store intersection points
		pointsVTKintersection = vtk.vtkPoints()
		code=obbTree.IntersectWithLine(pSource, pTarget, pointsVTKintersection, None)
		
		pointsVTKIntersectionData = pointsVTKintersection.GetData()
		noPointsVTKIntersection = pointsVTKIntersectionData.GetNumberOfTuples()
		pointsIntersection = []
		distanceCentroid=[]
		for idx in range(noPointsVTKIntersection):
			tup = pointsVTKIntersectionData.GetTuple3(idx)
			dx=(Cent[0]-tup[0])**2
			dy=(Cent[1]-tup[1])**2
			dz=(Cent[2]-tup[2])**2
			distanceCentroid.append((dx+dy+dz)**0.5)
			pointsIntersection.append(tup)	
			
		#Now find the furthest point from Centroid as Apex
		max_idx=np.argmax(distanceCentroid)
		apex_coord=pointsIntersection[max_idx]
		del obbTree,pointsVTKintersection,code,pointsVTKIntersectionData,noPointsVTKIntersection	
		return apex_coord

	def major_axis(self,Coords):
		pca=PCA(n_components=2)
		pca.fit(Coords)
		Norm1=pca.components_[0]
		Norm2=pca.components_[1]	
		return Norm1,Norm2

	def get_point_data(self,surface_data):
		#Number of points
		N_points=surface_data.GetNumberOfPoints()
		#Get Array Points
		Coords=np.zeros(shape=(N_points,3))
		Centroid=[0,0,0]
		for i in range(N_points):
			coord_=surface_data.GetPoint(i)
			Coords[i,0]=coord_[0]	
			Coords[i,1]=coord_[1]	
			Coords[i,2]=coord_[2]
			Centroid[0]+=coord_[0]
			Centroid[1]+=coord_[1]
			Centroid[2]+=coord_[2]
		Centroid[0]=Centroid[0]/N_points	
		Centroid[1]=Centroid[1]/N_points	
		Centroid[2]=Centroid[2]/N_points	
		return Centroid,Coords
