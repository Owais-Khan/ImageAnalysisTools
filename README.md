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
This script will convert a ".vti" image stack to a .raw file. This works seamlessly with the output of the ImageAnalysisLabelImage.py and ImageAnalysisLabelCenterlines.py (see above). You can chose to output a 0D array or a 3D array with the same shape as that of the original image stack. 
 
```console
foo@bar:~$ python ImageAnalysisVTKToRaw.py -InputFileName /path/to/input/image/file.vti -ReshapeArray 1
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



---


## 1.0 Dynamic CTA Perfusion Analysis
Dynamic CT Myocardial Perfusion Imaging (CT-MPI) is an advanced imaging modality that can be used to image the blood flow inside myocardial tissue. The image acquistion pipeline generally involves four steps: i) inducing vasodilator-mediate hypermia followed by injection of iodine-based contrast agent; ii) scanning myocardium for 20-30 seconds, providing ~15 time points during the passage of contrast; iii) extracting time attenuation curves in an artery (typically ascending aorta) and myocardium; and iv) calculating myocardial blood flow using a tracer kinetic model. 

![DynamicCTMPI Image1](images/1_DynamiCTMPI_Figure1.png)

*Figure 1: Stress dynamic CT-MPI shown at three of 15 time points to highlight the passage of contrast. The right panel shows arterial input function (AIF) and myocardial tissue time attenuation curves (TAC) that are needed in tracer kinetic model to compute myocardial blood flow. Green=Arterial input function, Red=Myocardial time attenuation curve.* 

### 1.1 Computing Coronary Centerlines
Coronary centerlines are needed to separate the myocardium into vessel-specific territories. The following spcript will take the "mesh-complete" folder from SimVascular as input and generate centerlines for each wall surface of coronary arteries. The centerlines will be stored inside "Centerlines" folder, the path to which you can assign as an argument. Please ensure to label the left and the right coronary trees as "wall_LCA\*.vtp" and "wall_RCA\*.vtp", respectively.

```console
foo@bar:~$ python ImageAnalysisCoronaryCenterlines.py -ifolder /path/to/mesh-complete/ -ofolder /path/to/outputfolder/to/store/centerlines/
```


### 1.2 Quantifying Vessel-Specific Myocardial Territories using Dynamic CT-MPI and coronary CTA
 You may use the following pipeline to separate the myocardium into territories based on proximity to a coronary vessel. These territorty maps are useful to determine vessel-specific ischemia, assign boundary conditions or compute vessel-specific myocardial mass. You will need the following three items to perform this analysis:
1. VTK Mesh that contains MBF data: This file should contain volumetric images of the myocardium and the scalar values of myocardial blood flow. The mesh should ideally be aligned/registered to the coronary CTA.
2. Centerlines folder that contains centerlines of all coronary arteries (i.e. the Centerlines folder you generated in Section 1.1)

```console
foo@bar:~$ python ImageAnalysisPerfusionCoronaryTerritories.py -ifile /path/to/MBF_Data.vtu/file -centerlinesfolder /path/to/Centerlines/folder -arrayname ImageScalars -ofile /path/to/MBF_Data_territories.vtu/file
``` 



 

