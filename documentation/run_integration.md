# Integration sub-routine

## Summary

This is the main function that executes sequentially all steps required to **integrate** (level 1a -> level 1b) a passive microwave radiometer. It has been designed for instruments using the hot/cold calibration scheme and primarly to GROMOS and SOMORA but has been adapted to MOPI 5 successfully and partly to MIAWARA-C which uses the tipping curve calibration scheme. 

Prior to the integration, a calibration step
([run_calibration](run_calibration.md)) must have been performed and saved for
this specific day and instrument otherwise this function will fail.

### Called from

[Main script](main.md)

### Calling

The following functions are called within the calibration functions. 

| name | inputs | outputs | type | Description |
|------|------|------|------|:-----------|
| read_level1a |calibrationTool |calibratedSpectra, meteoData, calibrationTool | Required | import calibrated data (must exist)
| add_meteo_data | calibrationTool, meteoData, calibratedSpectra | calibratedSpectra | Required | add meteo data to the calibrated spectra structure
| check_channel_quality | calibratedSpectra,calibrationTool,filterType | calibratedSpectra | Required | check the channel quality of the spectra within a structure array*
| window_correction | calibratedSpectra, calibrationTool | calibratedSpectra | Required | window correction for a spectrum*
| tropospheric_correction | calibratedSpectra, calibrationTool | calibratedSpectra | Required | tropospheric correction for a spectrum*
| integrate_calibrated_spectra | calibrationTool,calibratedSpectra | integratedSpectra | Required | integration of the calibrated spectra
| check_integrated | calibrationTool, integratedSpectra | integratedSpectra | Required | check of the integrated spectra and addition of some meta data
| plot_integrated_spectra | calibrationTool, integratedSpectra | - | Required | standard plot for level 1b
| save_level1b | calibrationTool, level1 | calibrationTool | Required | saves level 1b into netCDF file

\* these functions apply to both *calibratedSprectra* and *integratedSpectra*
structure array 

All functions are stored within the *calibrationTool* structure which is the single input for the *run_calibration* function.

## Input

*calibrationTool* with all required fields for performing an integration (see [*calibrationTool*](calibrationTool.md))

Also this function requires that a calibration was performed before for this day and the correspondings level 1a netCDF file.

## Outputs

Level 1b and more specifically:
* Standard plots in 1 PDF file
* Integrated spectra into a daily netCDF file (see [level 1](level1.md))

---
---

## Structure

To perform the integration, we executes the followings steps sequentially:

### 1. Importing the calibrated spectra (read_level1a)

These are the outputs from [*run_calibration*](run_calibration.md). Within the
GROSOM project, a generic function (read_level1_GROSOM.m) for reading netCDF file has been written
however, for speed purposes, we are using a quite dirty function reading only a
selection of parameters in the level 1a. Moreover, these parameters are still hardcoded in the code and this should be improved.

The imported calibrated spectra are then stored into a Matlab structure array, similarly as during the calibration routine (but with less variables). It is named *calibratedSpectra* and is then saved into a general *level1* structure which will eventually contain both the calibrated and integrated spectra at the end of the integration phase.

The meteo data contained in the level 1a are also read and stored within a standardized meteo structure as outputed by read_meteo_data during the calibration sub-routine.

### 2. Merge meteo data in *calibratedSpectra* (add_meteo_data)

As the meteo data are usually saved with their original time stamps into the level 1a, we need to merge them with the *calibratedSpectra* structure. For this, we use a generic function *add_meteo_data* which is simply integrating a classical *meteoData* structure into the *calibratedSpecta*. For most of the variables, this is a simple averaging on the defined *integrationTime*, except for the precipitation where we make a sum. 

---

### 3. Checking the spectrometer channels quality on the calibrated spectra (check_channel_quality)

In order to make some proper channel averages in the following functions, we need to identify the potential spurious channels contained at each time stamps of the *calibratedSpectra* structure array.

To do that, we have tested different kind of filters and we combine them with some empirical knowledge on the spectrometer channels (often, some channels are just always bad...). By specifying different filters to the *check_channel_quality*, it also enables to use a single function for checking the channel quality on the calibrated or on the integrated spectra: only a different filter is specified.

The different filter are selected by *filterType*, an integer that can take the following values:
1. a boxcar filter saved as *calibrationTool.filter1*
2. a boxcar filter saved as *calibrationTool.filter2*
3. a filter based on the standard deviation of individual channel to qualify the identify spurious channels. It uses *calibrationTool.maxStdDevTbCal* as a threshold and therefore is used for checks on *calibratedSpectra*.
4. a filter based on the standard deviation of individual channel to qualify the identify spurious channels. It uses *calibrationTool.maxStdDevTbInt* as a threshold and therefore is used for checks on *integratedSpectra*. 

As a general rule, we use filter n.3 on the calibrated spectra and a boxcar filter (n.2) to check the quality of the integrated spectra (see below)

### 4. Perform a window correction on calibrated spectra (window_correction)

The window correction accounts for the absorption of a potential window placed
between the instrument and the atmosphere (often the case for IAP radiometers at
least).

It computes the window (usually of styrofoam) brightness temperature with the
Planck law and considers that it absorb and re-emit part of the radiation coming
from outside at its own temperature.

In order to account for this effect, we need the transmission coefficient of the
window *calibrationTool.tWindow* (which should be re-computed regularly) as well
as its temperature (*TWindow*).

### 5. Perform a tropospheric correction on calibrated spectra (tropospheric_correction)

This function performs a tropospheric correction on the calibrated spectra. All
required parameters to define the correction are stored in
*calibrationTool.troposphericCorrection*, a sub-structure containing the
following the following parameters:
* *type*: type of tropospheric correction to apply ("Ingold_v1", "Ingold_v1_fit"
  or "Ingold_v2"). "Ingold_v1" and "Ingold_v2" use a mean tropospheric
  temperature computed from the *mean_air_temperature* (v1) or the *mean_air_temperature*
  and *mean_relative_humidity* (v2). "Ingold_v1_fit" uses 
* *useWings*: a string which specify which wing(s) are used to compute the mean
  wing brightness temperature ("both","left" or "right")
* *numberOfChannelsTropCorr*: number of channel to average together to compute a
  wing temperature.
* *skipFraction*: fraction of channel to keep at the end of the spectrum before
  computing the wing temperature.
* *deltaT*: constant temperature to remove rom the air_temperature to compute
  the mean tropospheric temperature (see Ingold et al.)

Note that if a window correction has been applied beforehand (which should be
the case), the tropospheric correction is made on the window corrected spectra.

At this stage, the goal here is not really to use the corrected spectra for the
integration but more to add some more knowledge (e.g. the opacity of the
atmosphere) before averaging the spectra together. In fact, depending on the
analysis that we want to do, this knowledge be used to select the calibrated
spectra that we integrate together.

### 6. Integrate the calibrated spectra (integrate_calibrated_spectra)

Function performing the integration of the calibrated spectra with different
options can be setup depending on the desired analysis. These options are listed
below and are set by two booleans into *calibrationTool*.
* Transmittance filtering of the calibrated spectra: only the spectra with an atmospheric transmittance higher
        than a threshold are kept for the integration.
* Flag filtering of the calibrated spectra: removing all flagged calibrated spectra before integration.

Once the required filtering of the spectra is done, a time averaging of the remaining
calibrated spectra contained within t_int is done. Note however that for some
variables, a sum is actually applied instead of a mean because more appropriate. This it
typically the case for the number of individual cycles or the amount of precipitation
recorded during t_int. The results of the integration are saved into a Matlab
structure array *integratedSpectra* within the global *level1* structure.

---

### 7. Checking the spectrometer channels quality on the integrated spectra (check_channel_quality)

Same as step 3 but applied on *integratedSpectra*.

### 8. Perform a window correction on integrated spectra (window_correction)

Same as step 4 but applied on *integratedSpectra*.

### 9. Perform a tropospheric correction on integrated spectra (tropospheric_correction)

Same as step 5 but applied on *integratedSpectra*.

### 10. Check the integrated spectra (check_integrated)

Similarly as during the calibration, this function checks the
*integratedSpectra* and adds some flags to the level 1b (see [quality
control](quality_control_calibration.md)). It also adds all required meta data
before the plotting and saving the data.

---

### 11. Plots integrated spectra (plot_integrated_spectra)

Standard plots for level 1b saved in PDF format. See the User Guide for a
detailed description of the plots. 

### 12. Save integrated spectra (save_level1b)

Saves the important variables into a single netCDF level 1b file (see [level
1](level1.md)). We keep a daily format and save all integration cycles of the
same day together. 

---
---

## Potential improvements

### Reading netCDF files in Matlab

Until now, quite ugly script for reading the level 1a during the integration
sub-routine. There is an automatic script existing (*read_level1_GROSOM.m*) but
it takes much longer to run (however, it reads all variables within level 1a/1b)

### Flag filtering of the calibrated spectra during integration 

Implement option to select which flags should be taken into account for the
integration.

### Check the window correction

Correction formula was taken from the old routines and might need a little check.