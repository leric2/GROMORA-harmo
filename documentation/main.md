# Main script

## Objective and role of this function

In the main script, we define all necessary parameters for the calibration and the integration.

While the main script is a good tool for research purposes, it will be adapted for an operational use.

### Called from

### Calling 
| name | type | Description |
|---|------|------|:-----------:|
| read_labview_log_generic | Optional |  | 
| import_default_calibration_tool | Required |  | 
| import_Instrument_calibrationTool | Required |  | 
| [run_calibration](run_calibration.md) | Required |  | 
| [run_integration](run_integration.md) | Required |  | 

## Inputs

| name | type | Description |
|---|------|------|:-----------|
| instrumentName | str | Name of the instrument to calibrate or integrate (case sensitive). 
| calibrationType | str | type of calibration to perform (standard or debug) 
| calibrate | boolean | defines if a calibration (level 0 to level 1a) is performed
| integrate | boolean | defines if an integration (level 1a to level 1b) is performed.
| readLabviewLog | boolean | defines if the labview log file should be read automatically
| labviewLogFolder | str | full path to locate the folder where is stored the labview log
| dates | datenum vector  | set of dates on which to perform the calibration.

Current options for the instrument names are:
* GROMOS
* SOMORA
* mopi5
* MIAWARA-C (to be checked)

Current status:

MIAWARA-C not working.

Note that for research purposes, some additional parameters can be modified,
either in the main script directly (see [Additional parameter](#define-some-additional-parameter-for-the-calibration)) or in
the import_Instrument_calibrationTool function.

## Structure

After setting the parameters, if a labview log file exist, the main script
begins by reading it and stores it into a *labviewLog* Matlab strucutre. This
will be further integrated within the *calibrationTool* structure.

After that, the main scripts begins to loop into the set of defined *dates* and
executes the following:

### 1. Import default *calibrationTool*

The main script is creates a generic *calibrationTool* structure containing
some common parameters for all instruments (mostly physical constants and time parameter for this day)

### 2. Define some additional parameter for the calibration

Mostly some boolean to decide for plotting or not some variables or some time related variables.

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

