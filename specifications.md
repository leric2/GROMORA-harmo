# Specifications
# DRAFT !
Collaborative and evolutive document defining all practical specifications for the SOMORA/GROMOS harmonized retrievals. It is focused on all computer and code related details for the harmonization project.

I would suggest to edit and approve a first version of this documents for the 1st part of this project (raw data -> calibrated spectra for both instruments, quality controlled) in 1-2 weeks from now (15.12.19). 

## Expectations and requirements
### Outputs
What we want as outputs for our retrievals:
* Calibrated spectra (level1)
* Ozone profile (level2)
* Full error characterization
* AVKs
* Quality flags easily checkables


### Sustainability
We want this code to be compatible with the following changes:
* Hardware
* ...
basic
### File format
Due to its widespread utilisation and its self-referencing property, we think that all level1 and level2 should be stored in the **netCDF** (binary) format. This option also enable to use all the existing tools to have a quick look at the data (Panoply,...)

The fact that the NDACC database requires hdf5 files is an additional point for this choice of data format as the 2 are almost similar. 

### Programming language 

The last available stable version of ARTS is used for this project. This is currently (22.11.19) version 2.3 which is 3 years old. 
The next ARTS stable version should be out Summer 2020 and therefore, we should try to have a code which is easily adaptable to a new version of ARTS (Jonas ?)

#### Working with ARTS:

Here are the 3 possible options to use ARTS for our retrievals (from Jonas presentation)
>
##### Use QPACK (MATLAB) → Input/Output are MATLAB structures
* High-Level interface (exact usage of ARTS is hidden)
* Provides OEM algorithms
* It’s what we have been using
>
##### Use the API (Python): Maps all WSMs and WSVs to Python functions / methods
* Directly calls ARTS
* High-Level interface(s) in developement
* It is the future.
>
##### Directly with .arts files (and therefore writing those files with either Matlab or Python) -> Input/Output are XML-files
* Very close to ARTS itself
* Manipulation of matrices is tedious
* Need code to write and parse XML data (e.g. to make plots)

The main programming language shall then be Matlab/Python because ...

## Timeline
Approximate timeline defined for this project (for now, focused on the first part of the project):

* 15.12.19: Specifications document edited and approved by everyone
* 31.12.19: Final draft ot the level0-level1 routines
* 15.01.20: Design of the quality control for level1 data.
* 15.02.20: Investigation of the level1 data for GROMOS and SOMORA
* 01.03.20: Start working on the level1 to level2 data
* ...


## Good practice
In order to have some coherence for this new retrieval code, here are some basic rules I would suggest to use when collaborating to the project.

### Modules
Every module has to make 1 action

### Variables 
No hardcoded variables, except in the main module


## Roles (self filling)

### Eliane

### Eric
Design and main contributor of the coding part
Organisation of the project
Git

### Gunter

### Klemens

### Axel



## Code main structure and data management  (maybe put that in another file)
**Here come some sketch for the main structure and data management that we are planning to use**

