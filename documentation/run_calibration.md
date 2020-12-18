# Calibration sub-routine

## Objective and role of this function

This is the main function that executes sequentially all steps required to calibrate (level 0 -> level 1a) a passive microwave radiometer. It has been designed for instruments using the hot/cold calibration scheme and primarly to GROMOS and SOMORA but has been adapted to MOPI 5 successfully and partly to MIAWARA-C which uses the tipping curve calibration scheme. 

### Called from

[Main script](main.md)

### Calling

The following functions are called within the calibration functions. 

| name | inputs | outputs | type | Description |
|------|------|------|------|:-----------|
| read_level0 | calibrationTool, readingRawFile | logFile, rawSpectra | Required | import raw files |
| harmonize_log | calibrationTool, logFile | logFile | Required | harmonize log structure
| reformat_spectra | rawSpectra, logFile, calibrationTool| rawSpectra | Conditional | format raw spectra into a matrix
| check_level0 | logFile, rawSpectra, calibrationTool | warningLevel0 | Optional | check raw files
| flip_spectra | rawSpectra | rawSpectra | Conditional | flip the raw spectra along the channels axis
| plot_raw_spectra | rawSpectra, plot parameters | - | Optional | basic plots of the raw counts
| read_meteo_data | calibrationTool | logFile.meteo | Required | read meteo data
| run_tipping_curve | rawSpectra, logFile, calibrationTool | logFile\. TC | Conditional | perform tipping curve calibration
| calibrate | rawSpectra, logFile, calibrationTool, calibrationTool.calType | drift, calibratedSpectra | Required | perform a hot-cold calibration
| check_calibrated | logFile, calibrationTool,calibratedSpectra | calibratedSpectra | Conditional | check calibrated spectra and add some metadata
| plot_calibrated_spectra | calibrationTool, drift, logFile.meteo, calibratedSpectra, N | - | Conditional | save standard plots for level 1a
| save_level1a | calibrationTool,logFile, calibratedSpectra, warningLevel0 | calibrationTool | Conditional | save level 1a in netCDF file

All functions are stored within the *calibrationTool* structure which is the single input for the *run_calibration* function.

## Inputs

*calibrationTool* with all required fields for performing a calibration (see [*calibrationTool*](calibrationTool.md))

## Outputs

Level 1a and more specifically:
* Plots in 2 PDF files
* Calibrated spectra into a daily netCDF file (see [level 1](level1.md))

## Structure

For the calibration of an instrument, we executes the followings steps sequentially:

### 1. Importing the raw files (read_level0)

Reading and formatting of the raw data for this day in 1 Matlab structure (log) and a Vector (binary).

This is a generic function unless there are multiple spectrometer data saved in
the same binary which makes it difficult to read in one go. 

All necessary parameters to locate and read the log and binary are stored in
*calibrationTool*. The *readingRawFile* file parameter can be use if one wants
to avoid reading the binary file for any reason.

### 2. Harmonization of the log file (harmonize_log)

While the log files are usually standardized text files, they have quite
different naming convention which, in addition, tend to vary with time for a
single instrument.

In order to keep the following of the routine generic, we therefore need to harmonize the
log file structure. 

This functions is therefore instrument specific and is basically only renaming
and combining some existing variables in the original log structure.

The standard outputs that are needed in the log structure (for a classical hot-cold calibration):

* time (datenum) and dateTime (datetime) variables constructed from the Year, Month, etc... 
* Year, Month, Day, Hour, Minute, Second
* Position
* Elevation_Angle
* Tipping_Curve_active, Tipping_Angle_Nr
* LN2_Sensors_OK, LN2_Level_OK: some flags extracted from (LN2_above_High, LN2_above_Low, LN2_Relay)
* T_Room
* T_Hot_Absorber (measured on the absorber)
* T_Window
* T_Out
* FE_T_Sys   	
* Ferranti_Lock	PLL_Lock	V_Gunn	
* FFT_adc_range, FFT_adc_overload, FFT_T_FPGA, FFT_Mode, FFT_Nr_of_acq	
* Data_file_size	
* SW_version, IWV (not required) 

If one or some variables do not exist in the original log, they are then created
in harmonize_log.

### 3. Format raw spectrum (reformat_spectra)
Transforms the rawSpectra from a vector to a matrix (if not already the case).

### 4. Raw data checks (check_level0)

Function to perfom some little checks on level0. This is also here that we
extract and save some extra timestamps that do not belong to this day.

The checks are quite basic and include:
1. Check that the number of bytes in the binary corresponds to the number of
   spectra recorded during the day and stored in the log file.
2. Check that the number of tipping curve makes sense
3. Check that the total number of individual cycle in the log file is realistic

After that, a check is made for detecting if all time stamps in the log file
actually belongs to the day we are calibrating. If this is not the case, it
saves the extra time stamps (log + spectra) into 2 extra raw files (see [quality_control_calibration](quality_control_calibration.md))

The results of this functions is a warning and a string describing the problems detected in
the file. This string is saved as an attributes in the level 1a.

### 5. Flip the spectrum (flip_spectra)

Some instrument actually records the spectra with an intermediate frequecy
[0:-500MHz and 0:500]. This functions just splits the whole raw matrix along the
the mid channel.

### 6. Plot the raw counts (plot_raw_spectra)

Simple function plotting (uglily) the raw counts. The three additional inputs
(lowerLim, upperLim, N) are just defining the lower/upper limit of the FFTS
counts and the number (N) of raw spectra to plot distributed regularly on the
whole day.

### 7. Import meteo data (read_meteo_data)

The second section of the calibration sub-routine is concerned with the meteorological
data needed for some instruments (typically using tipping curve calibration schemes). It
was decided to include it at this early stage of the calibration for all instruments
because the meteorological data can be useful to understand the radiometer observations.
We also decided that a few key meteorological parameters should be saved along the
calibrated spectra in the level 1a data files. 

Note that the meteo stations that are often located in the vicinity of the MWR
are usually different for every location. For this reason, there is an
instrument specific read_meteo_data function for each location. The output of
these functions are directly harmonized as a standard meteo structure containing
the following elements:
1. dateNum and dateTime: the time variables in Matlab datenum format (respectively datetime)
2. air_pressure: pressure of the air
3. air_temperature
4. rel_humidity: relative humidity at the station
5. precipitation: the rain accumulation measured between 2 time stamps in the original meteo file
6. tod: the time of day extracted from the meteo station time variable.

### 8. Run a tipping curve (run_tipping_curve)

Only for instruments working with tipping curve as main calibration scheme.

For more information, see the functions for MIAWARA-C written by Franzisca.

### 9. Perform a hot-cold calibration (calibrate)

In fact, the calibrate function basically applies the radiometric
calibration formula to the data within a given time interval
t_cal. It is a key function of the routine and where the actual calibration is made. 

See what is needed

### 10. Perform a hot-cold calibratin (check_calibrated)

I am then doing some other quality checks and adding flags in the function following the calibration before saving the calibrated spectra. 


### 11. Plots calibrated spectra (plot_calibrated_spectra)

Standard double page plots for level 1a. 

### 12. Save calibrated spectra (save_level1a)

Saves the important variables into a singel netCDF level 1a file (see [level 1](level1.md))
