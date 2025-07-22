import os
import numpy as np
import open3d as o3d
import argparse

class ImageAnalysisPointCloudToSurface():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFileName is None:
			self.Args.OutputFileName = os.path.splitext(self.Args.InputFileName)[0]
	def Main(self):
		print ("---Reading the Input File: %s"%self.Args.InputFileName)
		#Read the input data
		pcd = o3d.io.read_point_cloud(self.Args.InputFileName)

		#Estimate Normals
		print ("------ Estimating Normals")
		pcd.estimate_normals()

		#Using alpha shape function to get the geometry (lower alpha=smooth, higher alpha=more details)
		print ("------ Computing Geometry")
		mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, self.Args.Alpha)
		mesh.compute_vertex_normals()
		print ("------ Saving Surface: %s"%self.Args.OutputFileName+"_surface.ply")
		o3d.io.write_triangle_mesh(self.Args.OutputFileName+"_surface.ply", mesh)
	
		#Reading the DICOM Files
		

			
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take in a volumetric dataset with X Y Z T Mag U V W array written in matlab format.")
	parser.add_argument('-InputFileName', '--InputFileName', required=True, dest="InputFileName",help="File name contain coordinates in xyz format.")
	parser.add_argument('-Alpha', '--Alpha', required=False, dest="Alpha", default=2.0, help="Alpha value use to reconstruct surface from point cloud.")
	parser.add_argument('-OutputFileName', '--OutputFileName', required=False, dest="OutputFileName",help="File name to save the reconstructed surface.")
                        
	args=parser.parse_args()
	ImageAnalysisPointCloudToSurface(args).Main()

