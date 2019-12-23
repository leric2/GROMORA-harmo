# Calibration routine

## Objective and role of this function
Calibrate the observations from the instruments against the hot and cold load on a given "calibration time".


### Calls from
Main routine

### Calling
None

### Previously
Quality control of the raw data

### Next
Quality control for the calibrated spectra ? 

## Inputs
* Raw data
* Calibration time
* Instruments name ?

## Outputs
* Calibrated spectra
* Error code for the calibration

## Structure
Here, two main approaches can be taken:
1. An object-oriented approach:

Defining a set of tools that can are useful for calibrating a radiometer. For instance, we can think of defining different calibration method for the various instruments.

2. A more classical, step-by-step approach as was done previously

....

We should include the possibility to use the tipping curve calibrations done by the instruments approx. every 30 min (not exactly the same for GROMOS and SOMORA).


### Temperature

We have to make a decision regarding the Temperatures to be used for the calibration T_hot and T_cold. 

For T_hot, there is apparently 2 different measurements:
* AI_0 : in/near the absorber --> real temperature ?
* T_hot : near the heater, used for the stabilisation ?

For T_cold, we use the temperature of liquid nitrogen but we should also take into account some effect from the edge of the containter, reflections on the liquid surface, ...?. For these reason, T_cold was **80K** for GROMOS and only **77.15K** for SOMORA.

## How it was done in the past
### GROMOS


### SOMORA


## Questions remaining
