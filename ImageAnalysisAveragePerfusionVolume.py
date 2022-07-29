"""This script seperates perfusion data into individual 
images of each cycle and outputs the averaged image"""

import numpy as np
import glob
import os
import argparse
from vmtk import vmtkscripts

class CalculateAveragedImage():
    def __init__(self, Args):
        self.Args = Args
        filenames = glob.glob(f'{self.Args.InputFolderName}/*.dcm')
        filenames = sorted(filenames)
        self.N_file_per_cycle = int(len(filenames)/self.Args.NumberOfCycles)
        

        #self.N_img_per_cycle = int(len(filenames)/self.Args.NumberOfCycles)
        for i in range(0,self.Args.NumberOfCycles):
            directoryname = f'perfusion_image_cycle_{i}'
            pathDicom = f'{self.Args.InputFolderName}/{directoryname}'
            os.system(f"mkdir {pathDicom}")
            for j in range((i)*self.N_file_per_cycle,(i+1)*self.N_file_per_cycle-1):
                os.system(f'cp {filenames[j]} {pathDicom}')
            
            print(f'--- Looping over cycle: {i}')
            filenames_ = glob.glob(f'{pathDicom}/*.dcm')
            filenames_ = sorted(filenames_)
            os.system(f'vmtkimagereader -ifile {filenames_[0]} --pipe vmtkimagewriter -ofile {self.Args.InputFolderName}/Cycle_Image_{i}.vti')

    def main(self):
        Averaged_data = 0

        for i in range(0,self.Args.NumberOfCycles):
            
            image_name = f'Cycle_Image_{i}.vti'
            image_path = f'{self.Args.InputFolderName}/{image_name}'
            
            
            print(f'------------------------------------Cycle_Image_{i}')
                       
            reader = vmtkscripts.vmtkImageReader()
            reader.InputFileName = image_path
            reader.Execute()
            
            imageNumpyAdaptor = vmtkscripts.vmtkImageToNumpy()
            imageNumpyAdaptor.Image = reader.Image
            imageNumpyAdaptor.Execute()

            numpyImage = imageNumpyAdaptor.ArrayDict

            print(numpyImage)
            Averaged_data += np.array(numpyImage['PointData']['ImageScalars']) 
        
        Averaged_data = Averaged_data/13
        numpyImage['PointData']['ImageScalars'] = Averaged_data
        print(numpyImage)

        output_image = vmtkscripts.vmtkNumpyToImage()
        output_image.ArrayDict = numpyImage # ArrayDict_
        output_image.Execute()

        output_vti = vmtkscripts.vmtkImageWriter()
        output_vti.Image = output_image.Image
        output_vti.OutputFileName = f'{self.Args.InputFolderName}/{self.Args.OutputFileName}.vti'
        output_vti.Execute()

if __name__ == '__main__':
    #descreption
    parser = argparse.ArgumentParser(description='Thsi script takes a dicom folder with N cycles and outputs an averaged vti image')
    #Input
    parser.add_argument('-InputFolder', '--InputFolder', type = str, required = True, dest = 'InputFolderName', help = 'The name of the folder with all of the dicom files')
    #NumberOfCycles
    parser.add_argument('-NofCycle', '--NumberOfCycles', type = int, required = True, dest = 'NumberOfCycles', help = 'The number of perfusion images that are in the dicom folder')
    #OutputFileName
    parser.add_argument('-OutputFile', '-OutputFileName', type = str, required = True, default = 'Perfusionaveraged', dest = 'OutputFileName', help = 'The name of the output averaged file in vti format')
    args = parser.parse_args()
    CalculateAveragedImage(args).main()

    