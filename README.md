# GROMORA 

## Summary
Main repository for the harmonization of GROMOS and SOMORA ozone calibration and retrieval routines. To go in the direction of a unified calibration routine, the work from Franzisca on the MIAWARA-C instrument is now being incorporated to this repository as well. 

The goal of this repository is to enable anyone concerned at the IAP and at MCH to follow and comment on the progression of the new harmonized retrieval routines.

## Version
* 1.0 pre-release 12.2020
* 1.0 released 03.2021

## Structure of the repository
#### documentation 
Main folder containing the documentation of the whole routine, the specifications, the choices made for calibration, etc... It also contains a technical documentation for the functions controlling the calibration. 

#### scripts
Containing all scripts in 2 separated folders: 
1. calibration: level0 to level1b, in Matlab for now. 
2. retrieval: level1b to level2, in Python

Each folder contains additional specific documentation. 

#### sketch

Some sketch trying to explain the GROSOM routine on a graphical level

## Requirements

### calibration: 
Matlab 2020a or Matlab 2020b recommended. 

With a few adaptations (?), it should work with Matlab 2019 as well

### retrievals:
* ARTS 2.4
* Python 3.8
* pyretrievals

When the time comes, we will add some more information here.

## Contributions
Contributions from anyone are more than welcome according to the rules defined in the **contributions.md** documents.

