# ImageAnalysisTools
The scripts in this repository provide useful tools to analyze medical images. You will need to set up the following libraries for majority of these scripts:
1. vtk 
2. vmtk (see www.vmtk.org)
3. numpy
4. scipy

To get help for any script, please type:
```console
foo@bar:~$ python [ScriptName.py] -h
```
## Image Analysis Tools
The following scripts provide useful tools to perform Image Analysis operations:

### ImageAnalysisLabelImage.py
Label the Image (in .vti format) based on the Surface segmentation (in .vtp format). The surface must be closed with no holes or gaps. 

```console
foo@bar:~$ python ImageAnalysisLabelImage.py -InputFileName /path/to/volume/image.vti -InputSurface /path/to/segmented/surface.vtp -OutputFileName /path/to/outputimage.vti
``` 


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



 

