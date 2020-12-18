# Integration sub-routine

## Objective and role of this function

This is the main function that executes sequentially all steps required to calibrate (level 0 -> level 1a) a passive microwave radiometer. It has been designed for instruments using the hot/cold calibration scheme and primarly to GROMOS and SOMORA but has been adapted to MOPI 5 successfully and partly to MIAWARA-C which uses the tipping curve calibration scheme. 

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
* Plots in 1 PDF file
* Integrated spectra into a daily netCDF file (see [level 1](level1.md))

## Structure

During the integration, we executes the followings steps sequentially:

### 1. Importing the calibrated spectra (read_level1a)
saved from [*run_calibration*](run_calibration.md) 

### 2. Merge (add_meteo_data)



### 3. ... (check_channel_quality)

check the channel quality of the spectra within a structure array 

Both for calibrated and integrated spectra

HERE describe the different filters


### 4. ... (window_correction)

### 5. Perform a tropospheric correction (tropospheric_correction)

different option

### 6. Integrate the calibrated spectra (integrate_calibrated_spectra)

### 7. ... (check_channel_quality)

### 8. ...(window_correction)

### 9.Perform a tropospheric correction  (tropospheric_correction)

### 10. Perform a hot-cold calibratin (check_integrated)

### 11. Plots integrated spectra (plot_integrated_spectra)

### 12. Save integrated spectra (save_level1b)