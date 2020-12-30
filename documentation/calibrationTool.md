# calibrationTool structure

## Summary
This is the Matlab structure containing all the information to launch a calibration for
a MWR intrument. It contains not only all parameters required for calibrating a
specific instrument or day but also the functions that will be used to
perform the calibration.

---

## Table of Contents
1. [Building calibrationTool](#building-calibrationtool)
2. [Parameters](#parameters)
2. [Functions](#functions)

---

## Building calibrationTool

The *calibrationTool* structure is unique for each instrument and each date of processing. 

The building of the *calibrationTool* structure is done at different places
within the GROSOM routine:

### 1. *import_default_calibrationTool*

The creation of the *calibrationTool* is done within the *import_default_calibrationTool* function called from the [main script](main.md). After creation of the structure, *import_default_calibrationTool* also adds the variables that are not instrument specific, like the date and time variables and the physical constants. 

### 2. Main script

After *calibrationTool* has been initialized, some important parameters are defined directly within the [main script](main.md) of the GROSOM calibration. These are generally key parameters that are defined there because they are changing frequently. Most of them should be moved to the specific import functions when the routine will be ready for operationnal use.

### 3. *import_InstrumentName_calibrationTool*

This is where most of the fields from calibration Tool get filled. In this function, you have to setup all required parameters to perform the calibration and integration of this MWR. 

It is also there that are defined all the generic or specific functions that will be used in the [calibration](run_calibration.md) and [integration](run_integration.md) sub-routine. Remember that the both the calibration and integration routine only have *calibrationTool* as single input. 

At the end of *import_InstrumentName_calibrationTool*, there is also a place where we setup all the time dependent variables. 

### 4. During the execution

After the first 3 steps, the building of the *calibrationTool* structure is finished.

Only a few specific parameters will be added during the processing of the routine. 

---

## Parameters

### Time variables

|variable | type  | Description |
|------|------|:-----------|
| dateStr | str  | 'YYYY_MM_HH' |
| Year | double  | YYYY |
| Month | double  | MM |
| Day | double | DD |
| dateTime | datetime | Matlab datetime object for the day |  
| timeNumber | datenum | Matlab datenum of the day |  
| meanDatetimeUnit | str | the unit used for the time |  
| referenceTime | datenum | the time to take as reference for the level 1 time vector |  
| calendar | str | calendar type used |  
| timeZone | str | Matlab TimeZone |  
| calendar | str | type of calendar used for the time |  

---

### Metadata 
|variable | type | Description |
|-------|------|:-----------|
| instrumentName | str | name of the instrument |  
| dataLocation | str | location of the instrument | 
| PI_NAME | str | name of the PI | 
| PI_AFFILIATION | str | affiliation of the PI | 
| PI_ADDRESS | str | address of the PI | 
| PI_EMAIL | str | email of PI | 
| dataSource | str | MWR.O3_UBERN (see NDACC) | 


---


### Geolocation data 
|variable | type | unit | description | 
|------|------|------|:-----------|
| lon | double | degree_north | latitude defined according to WGS84
| lat | double | degree_east | longitude defined according to WGS84
| altitude | double | meter | altitude above see level
| azimuthAngle | double | degree measured clockwise positive, 0 deg is northwise | azimuth observation angle 


---


### Physical constants

| variable | type | unit | Description |
|------|------|------|:-----------|
| lightSpeed | double | m/s | speed of light |
| h | double | J*s | Planck constant |
| kb | double | J/K | Boltzmann constant  |
| zeroDegInKelvin | double | degree | conversion between degreeC and Kelvin |
| backgroundMWTb | double | K | MW background radiation |

---

### Spectrometer variables
| variable | type | Description |
|------|------|:-----------|
| observationFreq | double | observation frequency in \[Hz\] |  
| numberOfSpectrometer | double | number of spectro in this instrument |  
| spectro | str | name of the spectro |  
| samplingRateFFTS | double | temporal sampling rate of the spectro in \[MHz\] |  
| numberOfChannels | double | number of spectro channels |  
| IQProcessing | boolean | 1 for IQ processing spectro |  
| fLO1, fLO2, ... | double | value of the local oscillator frequencies in \[Hz\] |  
| LOFreqTot | double | total value of the local oscillator frequencies in \[Hz\] |  
| instrumentBandwidth | double | bandwidth of the spectro in \[Hz\] |  
| badChannels | double vector | containing individual channels of spectro known to work badly |  
| numberOfAquisitionSpectraHot | double | # of spectra aquired during a single hot measurement cycle |  
| numberOfAquisitionSpectraAntenna | double | # of spectra aquired during a single sky measurement cycle |  
| numberOfAquisitionSpectraCold | double | # of spectra aquired during a single cold measurement cycle  |  
| flipped_spectra | boolean | indicates if the spectra is flipped compared to normal |  
| flipAroundChannel | double | channel number where to split the spectra (usually mid-channel) |  

---

### Raw files check

| variable | type | Description |
|------|------|:-----------|
| checkLevel0 | boolean | 1 if we want to perform a check on the raw files | 
| tippingSize | double | expected size of the tipping curve (tc) cycle |  
| numberOfTippingCurveExpected | double | expected number of tc done in the day| 
| toleranceTippingCurves | double | tolerance for the # of tc during this day | 
| numberOfCyclesExpected | double | total number of cycle expected during this day | 
| toleranceNumberCycles | double | tolerance for the # cycle expected during this day  | 

---

### Files and folder variables 

Note: all folders needs to be given with their full paths.

| variable | type  | Description |
|------|------|:-----------|
| binaryDataExtension | str | extension of the binary raw file (usually  .bin) | 
| logFileDataExtension | str | extension of the binary log file (usually  .txt) | 
| bytesPerValue | double | number of bytes used by a single value in the binary files |  
| binaryType | str | ordering byte type of the binary file (often 'ieee-be') |  
| rawFileFolder | str | folder containing the raw and log files |  
| extraFileFolder | str | folder to use for storing extra time stamps |  
| level1Folder | str | folder to use for storing the outputs files (plots + level 1) |  
| meteoFolder | str | folder containing the meteo files | 
| filename | str | name of the raw file (often GROMOS09_dateStr) |  
| file | str |  rawFileFolder + filename|  
| delimiter_logfile | str  | delimiter symbol used in the log file (for instance '\t') |  
| labviewLog | boolean | 1 for trying to read the log book from the labview -> filename to specify |  
| filenameLevel1a | str | full name of level 1a netCDF output file |  
| filenameLevel1b | str | full name of level 1b netCDF output file |  


---

### Log variables 

These variables are used to understand the log file for this day. 

| variable | type  | Description |
|------|------|:-----------|
| THotUnit | str | unit of the recorded hot load temperature |  
| positionIndAsName | boolean | 1 if the position are stored as string (e.g. 'cold') instead of integer (0) in the log |  
| indiceCold | double | integer corresponding to position of the cold load |  
| indiceAntenna | double | integer corresponding to position of the sky observation |  
| indiceHot | double | integer corresponding to position of the hot load |  
| cycleDurationCold | double | approximate value of the cold measurement cycle time \[s\] | 
| cycleDurationSky | double | approximate value of the sky measurement cycle time \[s\] |  | 
| cycleDurationHot | double | approximate value of the hot measurement cycle time \[s\] |  | 

---

### Calibration variables, flags and outlier detection

The following variable are used for the calibration and the flagging of the
calibrated spectra. Most of the flagging variables and threshold are described
in more details in [quality control](quality_control_calibration.md).

| variable | type  | Description |
|------|------|:-----------|
| calType | str | type of the calibration to perform ('standard' or 'debug') | 
| outlierDectectionType | str | type of outlier detection to perform ('standard', 'noFFT' or 'none) |
| calibrationTime | double | time interval for the calibration \[min\] |  
| TCold | double | cold load temperature to use for the calibration |  
| TSysCenterTh | double | expected value of the noise receiver temperature |  
| TSysThresh | double | threshold value of the noise receiver temperature |  
| stdTSysThresh | double | threshold value of the standard deviation of the noise receiver temperature |  
| frequencyBandAroundCenterTSys | double | see [quality control](quality_control_calibration.md) |  
| THotTh | double | expected value of the hot load temperature |  
| THotAbsThresh | double | threshold value of the hot load temperature |  
| hotTemperatureStdThreshold | double | threshold value of the standard deviation of the hot load temperature |  
| stdAntAngleThresh | double  | threshold value of the standard deviation of the sky observation angle within a calibration cycle |  
| minNumberOfIndicePerCycle | double | minimum number of individual cycles to be averaged together for a valid calibration cycle |  
| maxProportionOfIndLN2LevelOutlier | double | see [quality control](quality_control_calibration.md) | 
| maxProportionOfIndLN2SensorOutlier | double | see [quality control](quality_control_calibration.md) | 



Used for outliers detection: 
| variable | type  | Description |
|------|------|:-----------|
| elevationAngleCold | double | expected elevation during cold measurement cycle |  
| elevationAngleAntenna | double | expected elevation during sky measurement cycle |  
| elevationAngleHot | double | expected elevation during hot measurement cycle |  
| elevationAngleColdTol | double | tolerance for outlier detection of single cold spectra based on elevation angle | 
| elevationAngleTolerance | double | tolerance for outlier detection of single sky spectra based on elevation angle | 
| elevationAngleHotTol | double | tolerance for outlier detection of single hot spectra based on elevation angle | 
| hotSpectraNumberOfStdDev | double | see [outlierDetection](#outlierdetection)|  
| coldSpectraNumberOfStdDev | double | see [outlierDetection](#outlierdetection) |  
| skySpectraNumberOfStdDev | double | see [outlierDetection](#outlierdetection) |  
| threshNumRawSpectraHot | double | see [outlierDetection](#outlierdetection) |  
| threshNumRawSpectraCold | double | see [outlierDetection](#outlierdetection) |  
| threshNumRawSpectraAnt | double | see [outlierDetection](#outlierdetection) | 

Note that all angles are in degree.

#### outlierDetection

During calibration process, one outlier detection technique consists at checking how consistent are the spectrometer counts on hot, cold and sky position within a day (or a calibration cycle in the case of sky observation). To do that, we check that every indidivual spectrum measured on the cold and hot position have a least a certain amount of channels (= threshNumRawSpectraHot, threshNumRawSpectraCold) that are comprised within:

Daily median +/- hotSpectraNumberOfStdDev * daily stdDev

In the case of the sky observation, we take into account a bigger variability of the atmosphere by using the median and stdDev of the calibration cycle (usually 10 minutes) to detect the potential outliers.

Note that in addition to this techniques, individual spectra can also be removed
by a check on the elevation angle or if the ADC of the spectro was overloaded
during recording. This is also explained in [quality control](quality_control_calibration.md).

---


### Meteo variables
| variable | type  | Description |
|------|------|:-----------:|
| doTippingCurve | boolean | 1 for instruments that performs a tipping curve calibration as main calibration scheme | 
| troposphericCorrection | structure | a structure containing all required information to perform a tropospheric correction (see [run_integration](run_integration.md)) | 

---

### Plotting variables

| variable | type  | Description |
|---|------|------|:-----------:|
| rawSpectraPlot | boolean | 1 to plot the raw counts (to be improved) | 
| calibratedSpectraPlot | boolean | 1 to plot the calibrated spectra | 
| calibratedSpectraSpectralPlot | boolean | 1 to plot the spectral variable for the calibrated spectra | 
| integratedSpectraPlot | boolean | 1 to plot the integrated spectra | 

---

### Integration variables


| variable | type  | Description |
|---|------|------|:-----------:|
| integrationTime | double | time interval for the integration \[min\] |  
| minNumberOfAvgSpectra | double | minimum number of calibrated spectra to get a valid integrated spectrum |  
| filterByTransmittance | boolean | decides if the calibrated spectra are filtered based on their atmospheric transmittance before integration |  
| transmittanceThreshold | double | transmittance threshold value for a "good" calibrated spectrum |  
| filterByFlags | boolean | decides if the calibrated spectra are filtered based on the level 1a flags before integration |  
| filterTypeChannelQualityCal | int | type of filtering use to identify spurious channels on the calibrated spectra |  
| filterTypeChannelQualityInt | int | type of filtering use to identify spurious channels on the integrated spectra |  
| filter1 | struct | see [run_integration](run_integration.md) |  
| filter2 | struct | see [run_integration](run_integration.md) |  
| maxStdDevTbCal | double | threshold for the standard deviation of individual channel on an calibration cycle | 
| maxStdDevTbInt | double | threshold for the standard deviation of individual channel on an integration cycle |  
| tWindow | double | microwave window transmittance value|  

For more information on these variables, see [run_integration](run_integration.md)

---

### Variable created during the calibration or integration

| variable | type  | Description |
|---|------|------|:-----------:|
| successfulCalibration | boolean | indicates a successful calibration |  
| flagVectorLength | double | number of level 1a flags (read from level 1a) |  
| logFile | struct | an sub structure where we stores some interesting meta data read from level 1a (like the flags meanings, comments, ...) |  
| successfulIntegration | boolean | indicates a successful integration |  

---
---

## Functions

Functions stored in calibrationTool only.

### Calibration function

All functions here are used in [run_calibration](run_calibration.md) and are detailed there.

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

---

### Integration function

All functions here are used in [run_integration](run_integration.md) and are detailed there.

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
