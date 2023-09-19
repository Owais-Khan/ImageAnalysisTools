The scripts in this repository are useful tools to analyze medical image data. You will need to set up the following libraries for majority of these scripts:
1. vtk 
2. vmtk (see www.vmtk.org)
3. numpy
4. scipy

To get help for any script, please type:
```console
foo@bar:~$ python [ScriptName.py] -h
```
# Image Analysis Tools
The following scripts provide useful tools to perform Image Analysis operations:

## ImageAnalysisDicomToVti.py
This script will take as an input a folder that contains a series of sub-folders, each with dicom files (for example, 4D Flow MRI, CT Perfusion). It will convert all of the image stacks in each of the sub-folder into a vti format file that can be read into paraview or analyzed using vtk library. You must also defined an output folder where are of the .vti files will be saved.
```console
foo@bar:~$  python ImageAnalysisDicomToVti.py -InputFolder /path/to/input/folder/ -OutputFolder ./path/to/output/folder/
```


## ImageAnalysisLabelImage.py
Label the Image (in .vti format) based on the Surface segmentation (in .vtp format). The surface must be closed with no holes or gaps. The output image will contain Labels=0 for outside surface and Labels=1 for inside surface.
```console
foo@bar:~$ python ImageAnalysisLabelImage.py -InputFileName /path/to/volume/image.vti -InputSurface /path/to/segmented/surface.vtp -OutputFileName /path/to/outputimage.vti


## ImageAnalysisLabelImage.py
Label the Image (in .vti format) based on the Surface segmentation (in .vtp format). The surface must be closed with no holes or gaps. The output image will contain Labels=0 for outside surface and Labels=1 for inside surface.
```console
foo@bar:~$ python ImageAnalysisLabelImage.py -InputFileName /path/to/volume/image.vti -InputSurface /path/to/segmented/surface.vtp -OutputFileName /path/to/outputimage.vti
``` 
## ImageAnalysisLabelCenterlines.py
Label the image (in .vti format) based on the centerline surface file (in .vtp format). The centerlines are obtained from VMTK. The output image will contain Labels=0 for all cells outside the centerline and Labels=1 for all cells on the centerline
```console
foo@bar:~$ python ImageAnalysisLabelCenterlines.py -InputFileName /path/to/volume/image.vti -InputSurface /path/to/centerlinesfile.vtp -OutputFileName /path/to/outputimage.vti
```
## ImageAnalysisVTKToRaw.py
This script will convert a ".vti" image stack to a .raw file. This works seamlessly with the output of the ImageAnalysisLabelImage.py and ImageAnalysisLabelCenterlines.py (see above). You can chose to output a 0D array or a 3D array with the same shape as that of the original image stack. You have to assigned the following options: -ArrayType, either as a PointData or CellData (use CellData for CenterlineLabels) and -ArrayName (use Labels for lumen and CenterlineLabels for centerline labels and -ArrayName (use Labels for lumen and CenterlineLabels for centerline labels)) 
 
```console
foo@bar:~$ python ImageAnalysisVTKToRaw.py -InputFileName /path/to/input/image/file.vti -ReshapeArray 1 -ArrayType CellData -ArrayName CenterlineLabels
```

## ImageAnalysisCenterlineFromSections.py
Compute a centerline by taking the centroid of each of the sections of a surface. We can control the smoothing of the centerline using XX parameters. The OutputFolder is optional. The centerline file will be stored in the OutputFolder/CenterlineFromSections.vtp for paraview-readable format and OutputFOlder/CenterlineFromSections.txt for an ascii format
```console
foo@bar:~$ python ImageAnalysisCenterlineFromSections.py -InputSurface /path/to/input/surface.vtp -OutputFolder /path/to/output/folder
```

## ImageAnalysisCenterlineDeviation.py
Compute the difference between two clipped and registered centerlines by finding the normal between the reference line and an input line. The deviation will be
stored in a paraview-readable format under CenterlineDeviation_Results/CenterlineDeviations.vtp. The OutputFolder is optional. 
```console
foo@bar:~$ python ImageAnalysisCenterlineDeviation.py -InputCenterline /path/to/input/centerline.vtp -BaseCenterline /path/to/reference/centerline.vtp -OutputFolder /path/to/output/folder
```
## ImageAnalysisProjectImageToMesh.py
This script can be used to project data from image/surface_mesh/volume_mesh onto another surface_mesh/volume_mesh. The -InputFileName1 is the source data and -InputFileName2 is the target mesh onto which the data must be projected. The -OutputFileName is optional and the script will use the addition of the two InputFileName as the OutputFileName.
```console
foo@bar:~$ python ImageAnalysisProjectImageToMesh.py -InputFileName1 /path/to/source/data/image/or/mesh -InputFileName2 /path/to/target/mesh/
```

## ImageAnalysisMovingAverageAlongCL.py
This script can be used to take an average of the scalar value along the centerline of a lumen.The scalar value should be stored as an array in the lumen volume meshfile. Alongside with the vtu file, the scripts need the centerline of the lumen and the surface of that. This script was originally written to extract the average pixel value along the centerline, which was the output of the ImageAnalysisProjectImageToMesh.py script. The Input arguments of the scripts are as follows: -InputVolumeFile takes the volumetric mesh including the scalar array, -InputCLFile takes the centerline of the lumen, -InputSurfaceFile takes the surface model file of the lumen. 
```console
foo@bar:~$ python ImageAnalysisMovingAverageAlongCL.py -InputVolumeFile /path/to/source/data/mesh -InputCLFile /path/to/target/centerline/file -InputSurfaceFile /path/to/target/model
```

## ImageAnalysisContrastDispersion.py
This script can be used to take the centerlines files of the lumen containing the averaged pixel value array and estimate the averaged blood velocity along the lumen. The centerline files are the outputs of Moving Average script. The analysis is better to be done on the straight part of the lumen and only the upslope samples in time domain. The Input arguments of the scripts is: -InputFolder takes the folder containing all of the CLFiles.
```console
foo@bar:~$ python ImageAnalysisContrastDispersion.py -InputFolder /path/to/folder/containing/clfiles 
```

## ImageAnalysisCreateStenosis.py
This scrip will create a stenosis in a given surface. You can pre-specify the location (e.g., 50% along the length of the centerline) and the diameter reduction (e.g., 0.5). The script takes a surface file and a centerlines file.
```console
foo@bar:~$ python ImageAnalysisCreateStenosis.py -InputFileName /path/to/input/surface -InputCenterlinesFileName -StenosisLocation 75
```
Optional Arguments:
- ```-StenosisDiameter```: Stenosis throat diameter as a fraction of the normal diameter (default=0.5).
- ```-StenosisLength```: Length of the stenosis as a function of the normal diameter (default=2).
- ```-OutputFileName```: Name of the output file (default is same as input file with stenosis tag).

---


## 1.0 Dynamic CTA Perfusion Analysis
Dynamic CT Myocardial Perfusion Imaging (CT-MPI) is an advanced imaging modality that can be used to image the blood flow inside myocardial tissue. The image acquistion pipeline generally involves four steps: i) inducing vasodilator-mediate hypermia followed by injection of iodine-based contrast agent; ii) scanning myocardium for 20-30 seconds, providing ~15 time points during the passage of contrast; iii) extracting time attenuation curves in an artery (typically ascending aorta) and myocardium; and iv) calculating myocardial blood flow using a tracer kinetic model. 

![DynamicCTMPI Image1](images/1_DynamiCTMPI_Figure1.png)

*Figure 1: Stress dynamic CT-MPI shown at three of 15 time points to highlight the passage of contrast. The right panel shows arterial input function (AIF) and myocardial tissue time attenuation curves (TAC) that are needed in tracer kinetic model to compute myocardial blood flow. Green=Arterial input function, Red=Myocardial time attenuation curve.* 


### 1.1 Projecting Volumetric MBF to Surface Maps.
MBF maps are volumetric and often challenging to visualize and interpret. We can project these volumetric maps on the surface of the myocardiam, which makes it easier to visualize the average MBF across each slice. The script will also fill any holes in the geometry with data interpolated from the neighbouring nodes. 

```console
foo@bar:~$ python ImageAnalysisMyocardiumProjectVolumeToSurface.py -InputVolume MBF_Volume.vtk -InputSurface Ventricle_Surface_CTA.vtp
```
Optional Arguments:
- ```-OutputSurface```: The filename used to store the output projected surface. If None, the output filename will be MBF_Volume_ProjectedMBF.vtp.

### 1.2 Compute Global Myocardial Blood Flow Statistics
You can output a text file that contains the following global statistics that can be used to normalize the MBF maps to reduce inter-patient variabilities.
- Average MBF
- Standard Deviation
- Statistical Mode
- 50th Percentile (Median)
- 75th Percentile
- 80% of 75th Percentile.
The script with output two files: i) ```MBF_Statistics.txt```, containing the above statistics, and ii) ```MBF_Raw.txt```, containing the raw values of MBF in an ASCII format. 

```console
foo@bar:~$ python ImageAnalysisMyocardiumStatistics.py -InputFileName /path/to/registered/MBF.vtk -OutputFolder /path/to/outputfolder/
```
Optional Arguments:
- ```-ArrayName```: The name of the array containing the MBF values in the vtk file. By default, it is ```scalars```.
### 1.3 Plotting Probability Density Function


### 1.4 Computing Coronary Centerlines
You can semi-automatically compute centerlines of all the vessles produced by SimVascular. These centerlines can be used to separate the myocardium into distinct coronary-specific territories. Note that the following script uses the VMTK package. When you run the script, it will loop through all the vessel surfaces in the input folder, and for each, ask you to define the Inlet and Outlet ids. Note: vessel wall surfaces must have the following L_*.vtp and R_*.vtp tags, which represents left and right coronary arteries.

```console
foo@bar:~$ python ImageAnalysisCoronaryCenterlines.py -InputFolderName /path/to/folder/with/wall/surfaces
```
Optional Arguments:
-```-OutputFolderName```: The output folder to store the centerlines. If None, the centerlines will be stored in /path/to/folder/with/wall/surfaces_Centerlines/.


### 1.5 Quantifying Vessel-Specific Myocardial Territories
You can separate the myocardium into territories based on coronary centerlines computed in Section 1.4. These territory maps can be used to compute vessel-specific ischemia, assign coronary boundary conditions and/or compute vessel-specific myocardial mass. You will need two input files:
- ```MBF.vtk```: Contains the MBF values in "scalar" array.
- ```CenterlinesFolder```: centerlines generated in Section 1.4
 
```console
foo@bar:~$ python ImageAnalysisMyocardiumTerritories.py -InputFileName /path/to/input/filename.vtk -CenterlinesFolder /path/to/Centerlines/folder -OutputFileName /path/to/MBF_Data_territories.vtu
``` 

Optional Arguments:
- ```-ArrayName```: Name of the array in the input file that contains MBF values. Default is "scalars".

### 1.6 Plotting Histogram of MBF Values
You can plot the histogram of the raw MBF values and compare it to normal distribution. This is one way to assess how skewed the probability distrubution funcition is compared to a gaussian normal distribution, and thus, assess the severity of ischemia. The script will also output the results of Shapri-Wilk normality test (with p<0.05 meaning it's not normally distributed).

```console
foo@bar:~$ python ImageAnalysisMyocardiumHistogramTecplot.py -InputFileName /path/to/input/filename.vtk -OutputFolder /path/to/output/folder/to/store/tecplotfile.dat
``` 

Optional Arguments:
- ```-ArrayName```: Name of the array in the input file that contains MBF values. Default is "scalars".
- ```-Bins```: Number of Bins for the histogram. Default is 300.
