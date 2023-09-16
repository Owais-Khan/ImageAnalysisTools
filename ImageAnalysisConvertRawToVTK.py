'''
Sep 13, 2023
Author: Anahita A. Seresti
The goal of this script is to take three 1D time-varying raw files and convert them 
into 3D vtk format files.
'''
import numpy as np
import vtk
import glob
import argparse
from utilities import numpy_to_vtk
class RawToVTK():
    def __init__(self,args):
        self.args = args
    def main(self):
        # Reading and Storing .raw datasets
        rawfile = glob.glob(f"{self.args.InputFolder}/*.raw")
        rawfile = sorted(rawfile)
        # Reading and Storing Header Information
        header = dict()
        with open(f"{self.args.InputFolder}/0.5_2.5_Descr.txt", "rt") as file:
            for line in file:
                header[f"{(line.split('<'))[1].split('>')[0]}"] = (line.split(">"))[1].split("<")[0]
        T = int(header["DimensionSizes TZYX"].split(" ")[0])
        dimZ = int(header["DimensionSizes TZYX"].split(" ")[1])
        dimY = int(header["DimensionSizes TZYX"].split(" ")[2])
        dimX = int(header["DimensionSizes TZYX"].split(" ")[3])
        VoxelSizeX = float(header["VoxelSizeX [mm]"])
        VoxelSizeY = float(header["VoxelSizeY [mm]"])
        VoxelSizeZ = float(header["VoxelSizeZ [mm]"])
        # Creating a vtkImageData to create an image with specified dimensions and voxel sizes
        data = vtk.vtkImageData()
        data.SetDimensions(dimX, dimY, dimZ)
        data.SetSpacing(VoxelSizeX, VoxelSizeY, VoxelSizeZ)
        data.SetOrigin(0, 0, 0)
        # Creating the vtk XML-type Image writer to save the image Information and execute the Image 
        writer = vtk.vtkXMLImageDataWriter()
        # Reading Velocity Data from raw files. The file has been written as float32 
        Xdata = np.fromfile(rawfile[0], dtype=np.float32)
        Ydata = np.fromfile(rawfile[1], dtype=np.float32)
        Zdata = np.fromfile(rawfile[2], dtype=np.float32)
        # Reshaping Velocity Data to its original shape
        Xdata = Xdata.reshape(T,dimZ,dimY,dimX)
        Ydata = Ydata.reshape(T,dimZ,dimY,dimX)
        Zdata = Zdata.reshape(T,dimZ,dimY,dimX)
        # Defining a vtk float array to store the velocity values
        vector = vtk.vtkFloatArray()
        vector.SetName('Velocity')
        vector.SetNumberOfComponents(3)
        vector.SetNumberOfTuples(data.GetNumberOfPoints())
        for t in range(0,T):
            #velocity_vector = np.empty((dimX,dimY,dimZ,3), dtype=np.float32)
            for x in range(dimX):
                for y in range(dimY):
                    for z in range(dimZ):
                        point_id = data.ComputePointId([x,y,z])
                        vector.SetTuple3(point_id, Xdata[t,z,y,x], Ydata[t,z,y,x], Zdata[t,z,y,x])
                #vector.append(np.array([Xdata[t,item], Ydata[t,item], Zdata[t,item]]))
            #vector = np.array(vector).reshape(dimX,dimY,dimZ,3)
            #print(np.array(vector).shape)
            print(f"--- Writing Image {t+1}")
            data.GetPointData().SetVectors(vector)
            writer.SetInputData(data)
            writer.SetFileName(f"{self.args.InputFolder}/Image{t}.vti")
            writer.Write()
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="The script takes a folder containing raw data and image header")
    parser.add_argument("-InputFolder", "--InputFolder", dest="InputFolder", type=str, required=True, help="Input the folder containing X,Y,Z raw dataset and header")
    args = parser.parse_args()
    RawToVTK(args).main()