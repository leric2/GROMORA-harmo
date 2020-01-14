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
* Liquid Nitrogen temperature to use ?

## Outputs
* Calibrated spectra in a structure

## Structure
For the calibration, I have started simple and decided to go with a step-by-step method similar to what was done in the past. For now, I implemented it as:
1. Selecting all the interesting cycles (c-a-h/h-a-c and no tipping curve) in a time interval (typically 10 minutes) and storing their start indices
2. Checking cold and hot spectra for wrong angles and spurious spectra. Spurious spectra are defined as a spectrum which has more than x channels that are not comprised in an interval of the median +/- 2*stdDev of this calibration cycle. 
3. We remove those spurious or wrong angle spectra and average the remaining hot and cold spectra for the time interval.
2. At this point, we also compute the mean hot temperature for this calibration cycle as well as its stdDev (flag over a certain threshold).
4. We then deals with the antenna indices where we check only for some bad pointing angles before doing the averaging and applying the calibration formula.
5. While we want only 1 brigthness temperature for this calibration cycle, we addtionnaly store cycle-by-cycle calibration as well as individual angle flags. 

This is what can be called the "calibration routine".

I am then doing some other quality checks and adding flags in the function following the calibration before saving the calibrated spectra. 

Everything is running for both instruments but the following points will have to be checked/discussed:

### To check:
* Are the calibrated spectra alright ? --> think about how to check thah now

### To discuss:
* For now, I have decided to use only **complete individual cycle** as valid cycle, is it sufficient ?
* What other quality flags do we want to add for the calibration cycle ?
* Tsys/std(Tsyss)/THot/std(THot)/numberOfAvgSpectra/... ?

We should also include the possibility to use the tipping curve calibrations done by the instruments approx. every 30 min (not exactly the same for GROMOS and SOMORA).

## Questions remaining
### Temperature

We have to make a decision regarding the Temperatures to be used for the calibration T_hot and T_cold. 

For T_hot, there is apparently 2 different measurements:
* AI_0 : in/near the absorber --> real temperature ?
* T_hot : near the heater, used for the stabilisation ?

There is also the following question: should we use T_hot only when the hot spectra are measured or should we use the average ot T_hot from all the measurements in the given calibration cycle. 

For T_cold, we use the temperature of liquid nitrogen but we should also take into account some effect from the edge of the containter, reflections on the liquid surface, ...?. For these reason, T_cold was **80K** for GROMOS and only **77.15K** for SOMORA.
