import numpy as np
import argparse
import skspatial.objects as SK


class ImageAnalysisCenterlinePlanarity():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFileName is None:
			self.Args.OutputFileName=self.Args.InputFileName.split(".")[-2]+"_Deviations.dat"

	def Main(self):
		#Read the Centerline 
		print ("--- Reading the Centerline Coordinates from: %s"%self.Args.InputFileName)
		Infile_=open(self.Args.InputFileName,'r')
		Points=[]
		Infile_.readline()
		for LINE_ in Infile_:
			line_=LINE_.split()
			Points.append([float(line_[0]),float(line_[1]),float(line_[2])])
		Infile_.close()
		Points=np.array(Points)
		print ("--- The Numbr of Centerline Points Read: %d"%len(Points))

		#Now compute the plane of best fit
		print ("--- Computing Plane of Best Fit")
		Plane = SK.Plane.best_fit(Points)
		print ("------ The Centroid of Plane is: ",Plane.point)
		print ("------ The Normal of Plane is:   ",Plane.normal)
		
		#Now compute the 
		print ("--- Computing the centerline distance from the plane of best fit")
		CL_Dev=[]
		CL_Dev_s=[]
		for Point_ in Points:
			CL_Dev.append(Plane.distance_point(Point_))
			CL_Dev_s.append(Plane.distance_point_signed(Point_))
		print ("--- Average Distance to Plane is: %.05f"%np.average(CL_Dev))
		
		#Now write the data to outputfile
		print ("--- Writing the data to output file: %s"%self.Args.OutputFileName)
		outfile=open(self.Args.OutputFileName,'w')
		outfile.write("Distance Distance_To_Plane Distance_To_Plane_signed\n")
		outfile.write("%.05f %.05f %.05f\n"%(0.0, CL_Dev[0],CL_Dev_s[0]))
		distance_=0
		for i in range(1,len(Points)):
			distance_+= np.linalg.norm(Points[i] - Points[i-1])
			outfile.write("%.05f %.05f %.05f\n"%(distance_,CL_Dev[i],CL_Dev_s[i]))	
		outfile.close()
		

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute the plane of best fit using centerline points and output a centerline deviation fill and average deviation. Note: Centerline in ascii format with 1st line as x,y,z header.")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the myocardial volume")
        
	#Output Filename 
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The filename to store the out of plane deviations.")

	args=parser.parse_args()
	ImageAnalysisCenterlinePlanarity(args).Main()
