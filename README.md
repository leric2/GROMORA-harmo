# GROMORA 

## Summary
Main repository for the harmonization of GROMOS and SOMORA ozone calibration and retrieval routines. To go in the direction of a unified calibration routine, the work from Franzisca and Alistair on the MIAWARA and MIAWARA-C instruments is now being incorporated to this repository as well. 

The goal of this repository is to enable anyone interested at the IAP and at MCH to follow, comment and contributes to the new harmonized retrieval routines.

## Version

### Calibration
* 1.0 pre-release 12.2020
* 1.0 released 03.2021
* 2.0 pre-released 06.2022
* 2.0 released 07.2022

### Retrievals:
* 1.0 never released
* 2.0 pre-released 06.2022
* 2.0 released 07.2022

## Structure of the repository
### documentation 
Main folder containing the documentation of the whole routine, the specifications, the choices made for calibration, etc... It also contains a technical documentation for the functions controlling the calibration and some sparse documentation of the retrievals routine. 

For more information on the routines, the user should read the calibration and retrievals user guides published on the BORIS database from the University of Berne:

* [Calibration routine for ground-based passive microwave radiometer: a user guide ](https://boris.unibe.ch/id/eprint/164418)
* [Harmonized ozone profile retrievals from GROMOS and SOMORA](https://boris.unibe.ch/170121/)

### scripts
Containing all scripts in separated folders: 
1. calibration: level0 to level1b, in Matlab for now. 
2. retrieval: level1b to level2, in Python
3. pyretrievals: the updated retrieval package, initially written by [Jonas Hagen](https://github.com/jonas-hagen/pyretrievals)
4. extra_scripts: some outdated script that might be useful once...

Each folder contains additional specific documentation. 

### sketch

Some sketch trying to explain the GROMORA routine on a graphical level (in progress)

## Requirements

### calibration: 
Matlab 2021a or Matlab 2021b recommended. 

With a few adaptations (?), it should work with Matlab 2019, 2020 and 2022 as well

### retrievals:
* ARTS 2.4 and the accompanying PyARTS package
* Python 3.8
* A list of required package listed in the [specifications](scripts/env_file_GROMORA.txt). 

## Operational retrievals

As all other instruments (at least until Summer 2023...), GROMOS operational retrievals are done on Stockhorn in the framework defined by Jonas (see [stockhorn-scripts](https://git.iap.unibe.ch/MW/stockhorn-scripts)). 

GROMOS provides both for radid delivery (RD) and the consolidated data to NDACC and therefore, each day is retrieved 2 times. The first time is done a few day after measurements using the operational ECMWF analysis data and a second time 6 month later using the ERA5 reanalysis. In principle, the RD data obtained after the first retrievals are sent directly to NDACC. After the second retrieval though, a human check should be performed before sending the consolidated data to NDACC.

As this will likely not happen in the future, this second retrieval and upload of the consolidated data has been automatize as well on Stockhorn. What happens is that there are two bash scripts running daily GROMOS retrievals and providing different retrieval date (few days ago or 6 months ago) and retrieval type (RD or consolidated) and both will upload directly the data to NDACC (once to the RD and once to the normal server).

Note that the calibration will be performed again after the 6 months, for the sake of space, we will only keep the latest calibrated data for each day.

## Results and data

The results of the first part of harmonization of the GROMORA time series have been submitted and will be linked here in due course.

The data from 2009 to 2022 for both instrument can be found on the [BORIS-portal](https://boris-portal.unibe.ch/cris/project/pj00023). Together with the data, you will find a full documentation of the resulting time series.

### Data analysis

There are 2 separate repositories which are used for the data analysis of the new harmonized GROMORA time series. 

* [Level1_Analysis](https://git.iap.unibe.ch/IAP_MCH/Level1_Analysis): not so nicely documented scripts to deal with the level 1, for instance to perform the concatenation of the level 1 files.
* [level2_Analysis](https://git.iap.unibe.ch/IAP_MCH/level2_analysis): the main repository for the analysis of the level 2 from GROMOS and SOMORA as well as all cross-comparisons with satellites.
