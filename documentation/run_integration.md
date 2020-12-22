# Integration sub-routine

## Objective and role of this function

This is the main function that executes sequentially all steps required to integrate (level 1a -> level 1b) a passive microwave radiometer. It has been designed for instruments using the hot/cold calibration scheme and primarly to GROMOS and SOMORA but has been adapted to MOPI 5 successfully and partly to MIAWARA-C which uses the tipping curve calibration scheme. 

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
| check_channel_quality | calibratedSpectra,calibrationTool,filterType | calibratedSpectra | Required | check the channel quality of the spectra within a structure array
| window_correction | calibratedSpectra, calibrationTool | calibratedSpectra | Required | window correction for a spectrum
| tropospheric_correction | calibratedSpectra, calibrationTool | calibratedSpectra | Required | tropospheric correction for a spectrum
| integrate_calibrated_spectra | calibrationTool,calibratedSpectra | integratedSpectra | Required | integration of the calibrated spectra
| check_integrated | calibrationTool, integratedSpectra | integratedSpectra | Required | check of the integrated spectra and addition of some meta data
| plot_integrated_spectra | calibrationTool, integratedSpectra | - | Required | standard plot for level 1b
| save_level1b | calibrationTool, level1 | calibrationTool | Required | saves level 1b into netCDF file

All functions are stored within the *calibrationTool* structure which is the single input for the *run_calibration* function.

## Inputs

*calibrationTool* with all required fields for performing an integration (see [*calibrationTool*](calibrationTool.md))

Also this function requires that a calibration was performed before for this day and the correspondings level 1a netCDF file.

## Outputs

Level 1b and more specifically:
* Standard plots in 1 PDF file
* Integrated spectra into a daily netCDF file (see [level 1](level1.md))

## Structure

To perform the integration, we executes the followings steps sequentially:

### 1. Importing the calibrated spectra (read_level1a)

These are the outputs from [*run_calibration*](run_calibration.md). Within the
GROSOM project, a generic function (read_level1_GROSOM.m) for reading netCDF file has been written
however, for speed purposes, we are using a quite dirty function reading only a
selection of parameters in the level 1a. Moreover, these parameters are still hardcoded in the code and this should be improved.

The imported calibrated spectra are then stored into a Matlab structure array, similarly as during the calibration routine (but with less variables). It is named *calibratedSpectra* and is then saved into a general *level1* structure which will eventually contain both the calibrated and integrated spectra at the end of the integration phase.

The meteo data contained in the level 1a are also read and stored within a standardized meteo structure as outputed by read_meteo_data during the calibration sub-routine.

### 2. Merge meteo data not *calibratedSpectra* (add_meteo_data)

As the meteo data are usually saved with their original time stamps into the level 1a, we need to merge them with the *calibratedSpectra* structure. For this, we use a generic function *add_meteo_data* which is simply integrating a classical *meteoData* structure into the *calibratedSpecta*. For most of the variables, this is a simple averaging on the defined *integrationTime*, except for the precipitation where we make a sum. 

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

### 5. Perform a tropospheric correction on calibrated spectra (tropospheric_correction)

different option

### 6. Integrate the calibrated spectra (integrate_calibrated_spectra)

### 7. Checking the spectrometer channels quality on the integrated spectra (check_channel_quality)

### 8. Perform a window correction on integrated spectra (window_correction)

### 9. Perform a tropospheric correction on integrated spectra (tropospheric_correction)

### 10. Check the integrated spectra (check_integrated)

### 11. Plots integrated spectra (plot_integrated_spectra)

### 12. Save integrated spectra (save_level1b)

## Potential improvements

### Reading netCDF files in Matlab

### Flag filtering of the calibrated spectra during integration 

Implement option to select which flags should be taken into account for the integration.