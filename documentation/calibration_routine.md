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
* Calibrated spectra in a structure
* Error code for the calibration

## Structure
For the calibration, I have started simple and decided to go with a step-by-step method similar to what was done in the past. For now, I implemented it sothat we can use two different methods:
* calibrate all individual cycles (2-1-0-0-1-2) of a given day
* regroup invidual cycle into "calibration cycle" every x minutes (for instance 10)

Both are running for both instruments but the following points will have to be checked/discussed:

### To check:
* Are the calibrated spectra alright ? --> think about how to check thah now

### To discuss:
* For now, I have decided to use only **complete individual cycle** as valid cycle, is it sufficient ?
* Angles are checked to include only spectra that have right pointing angle for hot/antenna/cold.
* Implementation of a quality check based on Tsys/std(Tsyss)/THot/std(THot)/numberOfAvgSpectra/... ?
* What do we need as outputs from the calibration ?
    * Tb
    * Time of the calibration
    * Number of individual (and valid) cycles considered in the calibration
    * Tsys and std(Tsys)
    * THot and std(THot)
    * Antenna pointing angle for every calibration/individual cycles ?
    * ...

We should also include the possibility to use the tipping curve calibrations done by the instruments approx. every 30 min (not exactly the same for GROMOS and SOMORA).

## Questions remaining
### Temperature

We have to make a decision regarding the Temperatures to be used for the calibration T_hot and T_cold. 

For T_hot, there is apparently 2 different measurements:
* AI_0 : in/near the absorber --> real temperature ?
* T_hot : near the heater, used for the stabilisation ?

For T_cold, we use the temperature of liquid nitrogen but we should also take into account some effect from the edge of the containter, reflections on the liquid surface, ...?. For these reason, T_cold was **80K** for GROMOS and only **77.15K** for SOMORA.
