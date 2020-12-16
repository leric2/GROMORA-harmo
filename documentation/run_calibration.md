# Calibration routine

## Objective and role of this function

### Calls from


### Calling


### Previously
...

### Next
...

## Inputs


## Outputs
* Calibrated spectra in a structure



## Structure

Reformat the spectrum

Optional

Transform the raw data (vector) into a matrix.

Level 0 checks

Optional

Perform a few checks on the raw data. Saved as an attributes in the level 1a.

\paragraph{Flip the spectrum}

Conditional

For some instrument, flipped spectrum

Inversion

\paragraph{Plot the raw counts}

Optional

Simple function plotting (uglily) the raw counts. The three additional inputs (lowerLim, upperLim, N) are just defining
the lower/upper limit of the FFTS counts and the number (N) of raw spectra to plot distributed regularly on the whole day.

## Outputs
* Calibrated spectra in a structure

## Structure
For the calibration, I have started simple and decided to go with a step-by-step method similar to what was done in the past. For now, I implemented it as:
1. Selecting all the indices (h-a-c, without tipping curve) in a time interval (typically 10 minutes)
2. Checking cold and hot spectra for wrong pointing angles and spurious spectra. Spurious spectra are defined as a spectrum which has more than x channels that are not comprised in an interval of the median +/- x*stdDev of this calibration cycle. 
3. We remove those spurious or wrong angle spectra and average the remaining hot and cold spectra for the time interval.
2. At this point, we also compute the mean hot temperature for this calibration cycle as well as its stdDev (flag over a certain threshold ?).
4. We then deals with the antenna indices where we check only for some bad pointing angles before doing the averaging and applying the calibration formula.
5. While we want only 1 brigthness temperature for this calibration cycle, we addtionnaly store cycle-by-cycle calibration as well as individual angle flags for debugging.

This is what can be called the "calibration routine".

I am then doing some other quality checks and adding flags in the function following the calibration before saving the calibrated spectra. 

### To check:
* Are the calibrated spectra alright ? --> think about how to check that now

### To discuss:
* What other quality flags do we want to add for the calibration cycle ?
* Tsys/std(Tsyss)/THot/std(THot)/numberOfAvgSpectra/... ?

## Questions remaining
### Temperature

We have to make a decision regarding the Temperatures to be used for the calibration T_hot and T_cold. 

For T_hot, there is apparently 2 different measurements:
* AI_0 : in/near the absorber --> real temperature ?
* T_hot : near the heater, used for the stabilisation ?

There is also the following question: should we use T_hot only when the hot spectra are measured or should we use the average ot T_hot from all the measurements in the given calibration cycle. 

For T_cold, we use the temperature of liquid nitrogen but we should also take into account some effect from the edge of the containter, reflections on the liquid surface, ...?. For