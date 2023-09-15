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
        rawfile = glob.glob(f"{self.args.InputFolder}/*.raw")
        rawfile = sorted(rawfile)
        header = dict()
        with open(f"{self.args.InputFolder}/1.5_2.5_Descr.txt", "rt") as file:
            for line in file:
                header[f"{(line.split('<'))[1].split('>')[0]}"] = (line.split(">"))[1].split("<")[0]
        T = int(header["DimensionSizes TZYX"].split(" ")[0])
        dimZ = int(header["DimensionSizes TZYX"].split(" ")[1])
        dimY = int(header["DimensionSizes TZYX"].split(" ")[2])
        dimX = int(header["DimensionSizes TZYX"].split(" ")[3])
        VoxelSizeX = float(header["VoxelSizeX [mm]"])
        VoxelSizeY = float(header["VoxelSizeY [mm]"])
        VoxelSizeZ = float(header["VoxelSizeZ [mm]"])
        data = vtk.vtkImageData()
        data.SetDimensions(dimX, dimY, dimZ)
        data.SetSpacing(VoxelSizeX, VoxelSizeY, VoxelSizeZ)
        data.SetOrigin(0, 0, 0)
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetInputData(data)
        Xdata = np.fromfile(rawfile[0], dtype=np.float32)
        Ydata = np.fromfile(rawfile[1], dtype=np.float32)
        Zdata = np.fromfile(rawfile[2], dtype=np.float32)
        Xdata = Xdata.reshape(dimX*dimY*dimZ,T)
        Ydata = Ydata.reshape(dimX*dimY*dimZ,T)
        Zdata = Zdata.reshape(dimX*dimY*dimZ,T)
        vector = []
        for t in range(0,T):
            for item in range(0, Xdata.shape[0]):
                vector.append([Xdata[item,t], Ydata[item,t], Zdata[item,t]])
                #vector = vector.reshape(dimX,dimY,dimZ,3)
            data.GetPointData().SetVectors(numpy_to_vtk(np.array(vector)))
            writer.SetFileName(f"{self.args.InputFolder}/Image{t}.vti")
            writer.Write()
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="The script takes a folder containing raw data and image header")
    parser.add_argument("-InputFolder", "--InputFolder", dest="InputFolder", type=str, required=True, help="Input the folder containing X,Y,Z raw dataset and header")
    args = parser.parse_args()
    RawToVTK(args).main()