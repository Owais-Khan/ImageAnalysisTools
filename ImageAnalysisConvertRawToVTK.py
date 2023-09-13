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
class RawToVTK():
    def __init__(self,args):
        self.args = args
    def main(self):
        textfile = glob.glob("{self.args.InputFolder}/*.txt")
        rawfile = glob.glob(f"{self.args.InputFolder}/*.raw")
        rawfile = sorted(rawfile)
        header = dict()
        with open(f"{self.args.InputFolder}/1.5_2.5_Descr.txt", "rt") as file:
            for line in file:
                header[f"{(line.split('<'))[1].split('>')[0]}"] = (line.split(">"))[1].split("<")[0]
        T = int(header["DimensionSizes TZYX"].split(" ")[0])
        Z = int(header["DimensionSizes TZYX"].split(" ")[1])
        Y = int(header["DimensionSizes TZYX"].split(" ")[2])
        X = int(header["DimensionSizes TZYX"].split(" ")[3])
        print(T, X, Y, Z)
        Xdata = np.fromfile(rawfile[0], dtype=np.float16)
        X_Vel = Xdata.reshape((X,Y,Z,T))
        Ydata = np.fromfile(rawfile[1], dtype=np.float16)
        Y_Vel = Ydata.reshape((X,Y,Z,T))
        Zdata = np.fromfile(rawfile[2], dtype=np.float16)
        Z_Vel = Zdata.reshape((X,Y,Z,T))
        print(X_Vel, Y_Vel, Z_Vel)
        
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="The script takes a folder containing raw data and image header")
    parser.add_argument("-InputFolder", "--InputFolder", dest="InputFolder", type=str, required=True, help="Input the folder containing X,Y,Z raw dataset and header")
    args = parser.parse_args()
    RawToVTK(args).main()