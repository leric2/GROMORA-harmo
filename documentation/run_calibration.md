# Calibration sub-routine

## Objective and role of this function

This is the main function that executes sequentially all steps required to **calibrate** (level 0 -> level 1a) a passive microwave radiometer. It has been designed for instruments using the hot/cold calibration scheme and primarly for GROMOS and SOMORA but has been adapted to MOPI 5 successfully and partly to MIAWARA-C which uses the tipping curve calibration scheme. 

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

All these functions are stored within the *calibrationTool* structure which is the single input for the *run_calibration* function.

## Input

*calibrationTool* with all required fields for performing a calibration (see [*calibrationTool*](calibrationTool.md))

## Outputs

Level 1a and more specifically:
* Plots in 2 PDF files
* Calibrated spectra into a daily netCDF file (see [level 1](level1.md))

---
---
## Structure

For the calibration of an instrument, the followings steps are executed sequentially:

### 1. Importing the raw files (read_level0)

Reading and formatting of the raw data for this day in 1 Matlab structure for the log text file and 1 Vector for the raw binary file.

This is a generic function unless there are multiple spectrometer data saved in
the same binary which makes it difficult to read in one go. 

All necessary parameters to locate and read the log and binary are stored in
*calibrationTool*. The *readingRawFile* file parameter can be use if one wants
to avoid reading the binary file for any particular reason.

### 2. Harmonization of the log file (harmonize_log)

While the log files are usually standardized text files, they have quite
different naming convention which, in addition, tend to vary with time for a
single instrument.

In order to keep the following of the routine generic, we therefore need to harmonize the
log file structure. 

This functions is therefore instrument specific and is basically only renaming
and combining some existing variables in the original log structure in order to get a standard *logFile* that can be used further.

The standard outputs that are needed in the log structure (for a classical hot-cold calibration):

* **time (datenum) and dateTime (datetime)** variables constructed from the Year, Month, etc...
* **Year, Month, Day, Hour, Minute, Second**
* **Position**
* **Elevation_Angle**
* **Tipping_Curve_active**, Tipping_Angle_Nr
* **LN2_Sensors_OK, LN2_Level_OK**: some flags extracted from (LN2_above_High, LN2_above_Low, LN2_Relay)
* T_Room
* **T_Hot_Absorber** (measured on the absorber)
* **T_Window**
* T_Out
* **FE_T_Sys**	
* Ferranti_Lock	PLL_Lock	V_Gunn	
* FFT_adc_range, **FFT_adc_overload**, FFT_T_FPGA, FFT_Mode, FFT_Nr_of_acq	
* Data_file_size	
* SW_version, IWV (not required) 

The required variables are in bold. If one or some of these variables do not
exist in the original log, they are then created in harmonize_log.

### 3. Format raw spectrum (reformat_spectra)
Transforms the rawSpectra from a vector to a matrix (if not already the case) with m line and n columns with:

* m: number of individual during this day
* n: number of channels on the spectrometer

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

Some instrument actually records the spectra with an intermediate frequency ranging from
[0:-500MHz and 0:500] (or an typical 1 GHz bandwidth). This function just mirrors the whole raw matrix along the
the mid channel.

### 6. Plot the raw counts (plot_raw_spectra)

Simple function plotting (uglily) the raw counts. The three additional inputs
(lowerLim, upperLim, N) are just defining the lower/upper limit of the FFTS
counts and the number (N) of raw spectra to plot distributed regularly on the
whole day.

--- 

### 7. Import meteo data (read_meteo_data)

The second section of the calibration sub-routine is concerned with the
meteorological data needed for some instruments already at the calibration level
(typically using tipping curve calibration schemes). It was decided to include
it at this early stage of the calibration for all instruments because the
meteorological data can be useful to understand the radiometric observations. We
also decided that a few key meteorological parameters should be saved along the
calibrated spectra in the level 1a data files so that it is easier to access it
if needed.

Note that the meteo stations that are often located in the vicinity of the MWR
are usually different for every location. For this reason, there is an
instrument specific read_meteo_data function for each location. The output of
these functions are directly harmonized as a standard meteo structure containing
the following elements:
1. *dateNum* and *dateTime*: the time variables in Matlab datenum (datetime) format 
2. *air_pressure*: absolute air pressure at the station
3. *air_temperature*
4. *rel_humidity*: relative humidity at the station
5. *precipitation*: the rain accumulation measured between 2 time stamps in the original meteo file
6. *tod*: the time of day extracted from the meteo station time variable.

### 8. Run a tipping curve (run_tipping_curve)

Only for instruments working with tipping curve as main calibration scheme.

For more information, see the functions for MIAWARA-C written by Franzisca.

---

### 9. Perform a hot-cold calibration (calibrate)

This is the main function of the calibration sub-routine. 

In fact, the calibrate function basically applies the radiometric
calibration formula to the data within a given time interval
t_cal. It is a key function of the routine and where the actual calibration is made. 

There are 2 calibration modes implemented within the main calibration function: "standard"
and "debug". In both modes, the calibration begins with a simplified channel averaged calibration which
is stored into a *drift* Matlab structure. It contains the following quantities for each
individual calibration cycle (h-a-c or c-a-h):


* Some time identifiers for this individual cycle
* Averaged raw counts on hot, cold and antenna positions
* Y-factor and receiver noise temperature 
* Calibrated brightness temperature (averaged on all channels)
* Daily median and standard deviation of hot, cold and sky spectra
* A list of outliers detected for this day (see [quality control](quality_control_calibration.md))


This first calibration output is used as a gross indication of data quality and
enable a good visual interpretation of the data for a given day. However, it
does not yet give a calibrated spectra computed on t_cal.

For the complete calibration, the following steps are executed sequentially:

1. Select observation indices that belong to this calibration cycle
2. Check individual cold and hot spectra for outliers based on pointing angles
   and other spurious spectra identification methods (see [quality control](quality_control_calibration.md))
3. Remove hot and cold outliers before taking the averaged hot and cold spectra for this calibration cycle
4. Compute T_hot, the spectral Y-factor and noise receiver temperature and
   additional parameter for this calibration cycle.
5. Identify spurious antenna spectra during this calibration cycle and remove them.
6. Use calibration equation to compute the sky brightness temperature using
   individual antenna spectrum and the mean hot and cold spectra.
7. Compute the mean brightness temperature and its standard deviation on $t_{cal}$.

Compared to the "standard" mode keeping only a single calibrated spectra for the
calibration cycle, the "debug" mode additionally performed individual cycle
spectral calibration. It also stores two averaged calibrated spectra computed
only during the h-a-c (respectively c-a-h) calibration scheme to see if the
order of target observation does impact the resulting calibration.

In both modes, the second output (after the *drift* structure) of the hot-cold
calibration function is a Matlab structure array with time interval t_{cal}
containing all the interesting parameters computed during the calibration
(Y-factor, T_hot, ...). Note that after this stage, we do not the raw data
matrix anymore.

### 10. check the calibration (check_calibrated)

In addition to the outliers removal during the calibration, there are some
further checks done on the calibrated data. These checks are done within this
function and result in the creation of a set of flags that are saved in the
level 1a (see [quality control](quality_control_calibration.md)). Contrary to
the outliers detection, the flags are indicative of the data quality and do not
lead to any data removal before level1a. Depending on the final use of the data,
the user is then free to take these flags into account of not. Also the flags
are determined for each calibration cycle (and not on individual spectrum).

When doing the checks on the calibrated spectra, we also add all the necessary
meta data added to the *calibratedSpectra* structure, extracted either from the
log file, either from *calibrationTool*. At this point, the calibration is then
completed we can begin the plotting and saving phase.

--- 

### 11. Plots calibrated spectra (plot_calibrated_spectra)

Standard page plots for level 1a saved in 2 PDF files. See the User Guide for a
detailed description of the plots. 

### 12. Save calibrated spectra (save_level1a)

Saves the important variables into a single netCDF level 1a file (see [level
1](level1.md)). We keep a daily format and save all calibration cycle of the
same day together. 
