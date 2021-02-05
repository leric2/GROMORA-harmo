# Main script

## Summary

In the main script, we define all required parameters for the calibration and
the integration of a given instrument. 

The main script is a tool for developement, especially to select some
dates and modify quickly some key variables. It also enable to use multiple
instrument which is not needed for an operationnal use. 

This is a good starting point to understand the functionning of the routine and
for running some specific dates.

### Calling 

Functions called directly from the main script: 

| name | type | Description |
|---|------|:-----------|
| read_labview_log_generic | Conditional | reads the labview text file (called once) 
| import_default_calibration_tool | Required | import the default *calibrationTool* structure 
| import_InstrumentName_calibrationTool | Required | build and complete the instrument specific *calibrationTool* structure
| [run_calibration](run_calibration.md) | Required | launch the calibration process 
| [run_integration](run_integration.md) | Required | launch the integration process 


## Parameters

The main parameters to launch the calibration and/or integration process within
the GROMORA project:

| name | type | Description |
|---|------|:-----------|
| instrumentName | str | Name of the instrument to calibrate or integrate (case sensitive). 
| calibrationType | str | type of calibration to perform (standard or debug) 
| calibrate | boolean | defines if a calibration (level 0 to level 1a) is performed
| integrate | boolean | defines if an integration (level 1a to level 1b) is performed.
| readLabviewLog | boolean | defines if the labview log file should be read automatically
| labviewLogFolder | str | full path to locate the folder where is stored the labview log
| dates | datenum vector  | set of dates on which to perform the calibration.

Current available options for the instrument names are:
* GROMOS: back to 10.03.2010 (meteo only back to 12.05.2010)
* SOMORA: in principle back to 2011 (at least) but in practice, lots of periods tried before 2015 have too much FFT overloads detected and deleted
* mopi5: 01.2019 - 06.2019
* MIAWARA-C: needs some additional work to implement Franzisca's work

Among the above mentioned, some periods might need some adjustements of the outlier detection to work properly.

Note that for research purposes, some additional parameters can be modified,
either in the main script directly (see the additional parameters below) or in
the *import_Instrument_calibrationTool function*.

## Structure

After setting the main parameters, if a labview log file exist, the main script
begins by reading it and stores it into a *labviewLog* Matlab structure. This
will be further integrated within the *calibrationTool* structure. Doing it first avoid reading it multiple times in the loop of dates.

After that, the main script begins to loop into the set of defined *dates* and
executes the following for each day:

### 1. Import default *calibrationTool*

The main script creates a generic *calibrationTool* structure containing
some common parameters (mostly physical constants and time parameter for this day which do not depend on the instrument)

### 2. Define some additional parameters for the calibration

Mostly some boolean to decide for plotting or not some variables or some time
related variables. In addition to these, we have kept some of the instrument
dependent parameters inside the main script to enable some quick access and
changes to a few key variables like the calibration and integration time (in
minutes) or the filtering options for the integration of the calibrated spectra
(see [calibrationTool](calibrationTool.md)).

### 3. Import instrument specific instrument parameters

Depending on *instrumentName*, some more parameters are defined which are
instrument specific. The rest of the parameters and functions to use for this
specific day and instrument are then read from a special import function which
is unique for each instrument. These are named
import\_instrumentName\_calibrationTool.m and contain all necessary parameters
for the rest of the processing. 

### 4. Execute calibration and integration

The last part of the main script is just calling the generic [run_calibration.m](run_calibration.md)
and [run_integration.m](run_integration.md) functions. Both functions only have a single input; the
*calibrationTool* structure and are the same for all instruments.

Note that if multiple spectrometers are present in a single instrument, we also
loop on the set of spectrometers for calibrating or integration all of them.

When the calibration and integration is over for this day, the next date starts.
Because the raw files are saved daily, we decided to keep daily reference.

## Some potential improvements

### Multi spectrometer capabilities 

Now saving specific level 1a and 1b files for each of the spectrometer. That is ok but could be even better if include all spectrometers data in 1 file at each level.

### Error messages

The transfer of error messages from the sub-routine *run_calibration* and *run_integration* can be improved.

### Operational script

Think about the operationnal implementation of this script including:
* files and folder path definition
* filtering and outliers removal
* ...

