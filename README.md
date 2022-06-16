# GROMORA 

## Summary
Main repository for the harmonization of GROMOS and SOMORA ozone calibration and retrieval routines. To go in the direction of a unified calibration routine, the work from Franzisca and Alistair on the MIAWARA and MIAWARA-C instruments is now being incorporated to this repository as well. 

The goal of this repository is to enable anyone interested at the IAP and at MCH to follow, comment and contributes to the new harmonized retrieval routines.

## Version

### Calibration
* 1.0 pre-release 12.2020
* 1.0 released 03.2021
* 2.0 released 06.2022

### Retrievals:
* 1.0 never released
* 2.0 released 06.2022

## Structure of the repository
#### documentation 
Main folder containing the documentation of the whole routine, the specifications, the choices made for calibration, etc... It also contains a technical documentation for the functions controlling the calibration and some sparse documentation of the retrievals routine. 

For more information on the routines, the user should read the calibration and retrievals user guides located in a separated repository [UserGuideGROMORA](https://git.iap.unibe.ch/IAP_MCH/UserGuideGROMORA.git) or published on the BORIS database from the University of Berne.

#### scripts
Containing all scripts in 2 separated folders: 
1. calibration: level0 to level1b, in Matlab for now. 
2. retrieval: level1b to level2, in Python

Each folder contains additional specific documentation. 

#### sketch

Some sketch trying to explain the GROMORA routine on a graphical level (in progress)

## Requirements

### calibration: 
Matlab 2020a or Matlab 2020b recommended. 

With a few adaptations (?), it should work with Matlab 2019 as well

### retrievals:
* Linux
* ARTS 2.4 and the accompanying PyARTS package
* Python 3.8
* A list of required package listed in the [specifications](scripts/env_file_GROMORA.txt). 

## Results and data

The results of the first part of harmonization of the GROMORA time series have been submitted and will be linked here in due course.

The data from 2009 to 2022 for both instrument can be found on the BORIS-portal (LINK HERE). Along with the data, you will find a full documentation of the resulting time series.

### Data analysis

There are 2 separate repositories which are used for the data analysis of the new harmonized GROMORA time series. 

* [Level1_Analysis](https://git.iap.unibe.ch/IAP_MCH/Level1_Analysis): not so nicely documented scripts to deal with the level 1, for instance to perform the concatenation of the level 1 files.
* [level2_Analysis](https://git.iap.unibe.ch/IAP_MCH/level2_analysis): the main repository for the analysis of the level 2 from GROMOS and SOMORA as well as all cross-comparisons with satellites.

## Contributions
Contributions from anyone are more than welcome according to the rules defined in the **contributions.md** documents.

