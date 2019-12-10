
# DRAFT
Collaborative and evolutive document defining all practical specifications for the SOMORA/GROMOS harmonized retrievals. It is focused on all computer and code related details for the harmonization project.

I would suggest to edit and approve a first version of this documents for the 1st part of this project (raw data -> calibrated spectra for both instruments, quality controlled) in 1-2 weeks from now (15.12.19).

## Expectations and requirements
### File format
Due to its widespread utilisation and its self-documentation property, we think that all level1 and level2 should be stored in the **netCDF** (binary) format. This option also enable to use all the existing tools to have a quick look at the data (Panoply,...)

The fact that the NDACC database requires hdf5 files is an additional point for this choice of data format as the 2 are almost similar. 

### Level of data
#### Level 0
Those are the raw data and for both GROMOS and SOMORA, are composed of:
1. *.csv* file containing all the meta data for each measurements
2. *.bin* file containing the raw measured spectra

#### Level 1a: calibrated spectra
For every calibration cycle, we store and save the raw calibrated spectra in a dedicated netCDF file (without removing outliers). Also, according to the stability of the hot and cold load, we might consider defining a "accumulation" or "pre-integration" time for calibrating the specta (which then has to be noted in the meta data).

#### Level 1b: integrated and corrected calibrated spectra
Depending on the anaylsis, we might need different integration time and take different decision regarding what to do with the outliers. It might make sense then to do this now and not before storing the calibrated spectra. In this step, we would also include all correction that are instrument dependent.

#### Level 2: Ozone profile
The final product for our retrieval. It will also be stored as netCDF file and might include some other retrieved quantities like water-vapor for instance (if we decide to include the tropospheric correction in the retrieval).

### Additional outputs
What we want as additionnal outputs for our retrievals:
* Full error characterization
* AVKs
* Quality flags for every levels of data
* Other products ?

### Sustainability and reproductibility
We want this code to be compatible with changes in:
* Hardware<br/> 
Especially, we have to keep in mind the changes in spectrometer to include easily the older data from the instruments.

* Software<br/>
New version of Matlab<br/>
New version of ARTS 

### Programming language 
The main programming language shall then be **Matlab** because this is the language currently in use at the IAP. 

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

## Timeline
Approximate timeline defined for this project (for now, focused on the first part of the project):

* 15.12.19: Specifications document edited and approved by everyone
* 31.12.19: Final draft ot the level0-level1a routines
* 15.01.20: Design of the quality control for level1 data.
* 15.02.20: Investigation of the level1 data for GROMOS and SOMORA
* 01.03.20: Start working on the level1 to level2 data
* ...

## Good practice
In order to have some coherence for this new retrieval code, here are some basic rules I would suggest to use when collaborating to the project.

### Modules
Every module has to make 1 (and only 1) action at a time and be named accordingly. Suggestions for naming the modules:

*module_doing_this_action(var1,var2,...)*

Another question to be solved regarding the modules is:

*Do we want to have one specific module for each instruments when needed (corrections, etc...) or do we prefer to have some switch case and the instrument's name as input ?*

### Variables 
No hardcoded variables, except in the main module. Suggestions for naming the variables:

*nameOfTheVariable*

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
![MockUp for the project] (https://git.iap.unibe.ch/IAP_MCH/GROSOM-harmo/src/branch/master/sketch/MockUp_GROSOM.png)