# Specifications
Collaborative document defining all practical specifications for the SOMORA/GROMOS harmonized retrievals. It is focused on all computer and code related details for the harmonization project.

I would suggest to edit and approve a first version of this documents for the 1st part of this project (raw data -> calibrated spectra for both instruments, quality controlled). 

## Expectations and requirements
### Outputs
What we want as outputs for our retrievals:
* Calibrated spectra (level1)
* Ozone profile (level2)
* Full error characterization
* AVK
* Quality flags



### Sustainability
We want this code to be compatible with the following changes:
* Hardaware
* ...

### File format
Due to its widespread utilisation and its self-referencing property, we think that all level1 and level2 should be stored in the **netCDF** (binary) format. This option also enable to use all the existing tools to have a quick look at the data (Panoply,...)

The fact that the NDACC database requires hdf5 files is an additional point for this choice of data format as the 2 are almost similar. 

### Programming language 

The last available stable version of ARTS is used for this project. This is currently (22.11.19) version 2.3 which is 3 years old. 
The next ARTS stable version should be out Summer 2020 and therefore, we should try to have a code which is easily adaptable to a new version of ARTS (Jonas ?)

#### Working with ARTS:

In my opinion, there are 3 possible ways to use ARTS for our retrievals:
* Matlab with QPack
* Python with pyretrievals and Typhon package
* Directly with .arts files


## Timeline
Approximate timeline defined for this project (for now, focused on the first part of the project):

* 15.12.19: Specifications document edited and approved by everyone
* 31.12.19: Final draft ot the level0-level1 routines
* 15.01.20: Design of the quality control for level1 data.
* 15.02.20: Investigation of the level1 data for GROMOS and SOMORA







## Good practice
In order to have some "coding" coherence for this new retrieval, here are some general 
### Modules
Every modules

### Variables 

## Roles

### Eliane

### Eric
Design and main contributor to the coding part. 



## Code main structure and data management
