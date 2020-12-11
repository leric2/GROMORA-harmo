# calibrationTool structure

## Summary
This is the structure containing all the information to launch a calibration of GROSOM intrument. It contains not only all parameters required for calibrating a given instrument or a given day
but also the functions that will be used to performed the calibration.

All parameters must contain: type, units, where it is defined


## Parameters

### Time variable 

| # | variable | type  | Description |
|---|------|------|:-----------:|
| 1 | dateStr | str  | 'YYYY_MM_HH' |
| 2 | Year | double  | YYYY |
| 3 | Month | double  | MM |
| 4 | Day | double | DD |
| 5 | dateTime | datetime |  |  
| 6 | timeNumber | datenum |  |  
| 7 | meanDatetimeUnit |  |  |  
| 8 | referenceTime |  |  |  
| 9 | calendar |  |  |  
| 10 | timeZone |  |  |  
| 11 | calendar |  |  |  

### Metadata 

| # | variable | type | Description |
|---|------|------|:-----------:|
| ... | instrumentName | str |  |  
| 11 | dataLocation |  |  | 
| 11 | PI_NAME |  |  | 
| 11 | PI_AFFILIATION |  |  | 
| 11 | PI_ADDRESS |  |  | 
| 11 | PI_EMAIL |  |  | 
| 11 | dataSource |  |  | 


### Geolocation data
 | # | variable | type | Description |
|---|------|------|:-----------:|
| 11 | lon |  |  | 
| 11 | lat |  |  | 
| 11 | altitude |  |  | 
| 11 | azimuthAngle |  |  | 

### Physical constant / parameters

| # | variable | type | unit | Description |
|---|------|------|------|:-----------:|
| 8 | lightSpeed | double | m/s |  |
| 9 | h | double |  |  |
| 10 | kb | double |  |  |
| 11 | zeroDegInKelvin | double |  |  |
| 12 | backgroundMWTb | double |  | MW background radiation |
| 12 | backgroundMWTb | double |  | MW background radiation |

### Spectrometer variables
| # | variable | type | Description |
|---|------|------|:-----------:|
| ... | observationFreq | double |  |  
| ... | numberOfSpectrometer | double |  |  
| ... | spectrometer |  |  |  
| ... | samplingRateFFTS | double |  |  
| ... | numberOfChannels | double |  |  
| ... | IQProcessing |  |  |  
| ... | fLO1, fLO2, ... | double |  |  
| ... | LOFreqTot | double |  |  
| ... | instrumentBandwidth | double | bandwidth \[Hz\] |  
| ... | badChannels |  |  |  
| ... | numberOfAquisitionSpectraHot |  |  |  
| ... | numberOfAquisitionSpectraAntenna |  |  |  
| ... | numberOfAquisitionSpectraCold |  |  |  
| ... | flipped_spectra | boolean  |  |  

### Raw files check

| # | variable | type | Description |
|---|------|------|:-----------:|
| ... | checkLevel0 | boolean |  | 
| ... | tippingSize | double |  |  
| ... | numberOfTippingCurveExpected | double |  | 
| ... | toleranceTippingCurves | double |  | 
| ... | numberOfCyclesExpected | double |  | 
| ... | toleranceNumberCycles | double |  | 

### Files and folder variables
| # | variable | type  | Description |
|---|------|------|:-----------:|
| ... | binaryDataExtension | str |  | 
| ... | logFileDataExtension |  |  |  
| ... | bytesPerValue |  |  |  
| ... | binaryType |  |  |  
| ... | rawFileFolder |  |  |  
| ... | extraFileFolder |  |  |  
| ... | level1Folder |  |  |  
| ... | meteoFolder | str |  | 
| ... | filename |  |  |  
| ... | file |  |  |  
| ... | delimiter_logfile |  |  |  
| ... | labviewLog | boolean |  |  
| ... | filenameLevel1a |  |  |  
| ... | filenameLevel1b |  |  |  


### Log variables

| # | variable | type  | Description |
|---|------|------|:-----------:|
| ... | THotUnit | str |  |  
| ... | positionIndAsName | boolean |  |  
| ... | indiceCold | double |  |  
| ... | indiceAntenna | double |  |  
| ... | indiceHot | double |  |  
| ... | elevationAngleCold | double |  |  
| ... | elevationAngleAntenna | double |  |  
| ... | elevationAngleHot | double |  |  
| ... | elevationAngleColdTol | double |  | 
| ... | elevationAngleTolerance | double |  | 
| ... | elevationAngleHotTol | double |  | 
| ... | cycleDurationCold | double |  | 
| ... | cycleDurationSky | double |  | 
| ... | cycleDurationHot | double |  | 

### Calibration variables + outlier detection ?
| # | variable | type  | Description |
|---|------|------|:-----------:|
| ... | calType | str |  | 
| ... | calibrationTime |  | s |  
| ... | TSysCenterTh |  |  |  
| ... | TSysThresh |  |  |  
| ... | stdTSysThresh |  |  |  
| ... | frequencyBandAroundCenterTSys |  |  |  
| ... | THotTh |  |  |  
| ... | THotAbsThresh |  |  |  
| ... | hotTemperatureStdThreshold |  |  |  
| ... | stdAntAngleThresh |  |  |  
| ... | minNumberOfIndicePerCycle |  |  |  
| ... | hotSpectraNumberOfStdDev |  |  |  
| ... | coldSpectraNumberOfStdDev |  |  |  
| ... | threshNumRawSpectraHot |  |  |  
| ... | threshNumRawSpectraCold |  |  |  
| ... | threshNumRawSpectraAnt |  |  | 
| ... | maxProportionOfIndLN2LevelOutlier |  |  | 
| ... | maxProportionOfIndLN2SensorOutlier |  |  | 
| ... | threshNumRawSpectraAnt |  |  |
| ... | maxStdDevTbCal |  |  |  
| ... | goodFlagLN2Above |  |  |  
| ... | goodFlagLN2Below |  |  |  

### Meteo variables
| # | variable | type  | Description |
|---|------|------|:-----------:|

* doTippingCurve
* tWindow
* troposphericCorrection

### Plot variables

| # | variable | type  | Description |
|---|------|------|:-----------:|
| ... | rawSpectraPlot | boolean |  | 
| ... | calibratedSpectraPlot | boolean |  | 
| ... | calibratedSpectraSpectralPlot | boolean |  | 
| ... | integratedSpectraPlot | boolean |  | 

### Integration variables
| # | variable | type  | Description |
|---|------|------|:-----------:|
| ... | integrationTime |  | s |  
| ... | minNumberOfAvgSpectra | double |  |  
| ... | filterByTransmittance | boolean |  |  
| ... | filterByFlags | boolean |  |  
| ... | filterTypeChannelQualityCal | int |  |  
| ... | filterTypeChannelQualityInt | int |  |  
| ... | filter1 | struct |  |  
| ... | filter2 | struct |  |  
| ... | maxStdDevTbInt | boolean |  | 

### Variable list
* successfulCalibration
* flagVectorLength
* logFile
* 
* successfulIntegration

## Functions

### Main scripts function
| # | name | type | Description |
|---|------|------|:-----------:|
| 1 | read_labview_log_generic | Optional |  | 
| 2 | import_default_calibration_tool | Required |  | 
| 3 | import_Instrument_calibrationTool | Required |  | 
| 4 | run_calibration | Required |  | 
| 5 | run_integration | Required |  | 


### Calibration function
| # | name | inputs | outputs | type | Description |
|---|------|------|------|------|:-----------:|
| 6 | read_level0 | calibrationTool, readingRawFile | logFile, rawSpectra | Required | import raw files |
| 7 | harmonize_log | calibrationTool, logFile | logFile | Required | ... |
| 8 | reformat_spectra | rawSpectra, logFile, calibrationTool| rawSpectra | Conditional | ... |
| 9 | check_level0 | logFile, rawSpectra, calibrationTool | warningLevel0 | Optional | ... |
| 10 | flip_spectra | rawSpectra | rawSpectra | Conditional | ... |
| 11 | plot_raw_spectra | rawSpectra, plot parameters | - | Optional | ... |
| 12 | read_meteo_data | calibrationTool | logFile.meteo | Optional ?| ... |
| 13 | run_tipping_curve | rawSpectra, logFile, calibrationTool | logFile.TC | Conditional | ... |
| 14 | calibrate | rawSpectra, logFile, calibrationTool, calibrationTool.calType | drift, calibratedSpectra | Required | ... |
| 15 | check_calibrated | logFile, calibrationTool,calibratedSpectra | calibratedSpectra | Conditional | ... |
| 16 | plot_calibrated_spectra | calibrationTool, drift, logFile.meteo, calibratedSpectra, N | - | Conditional | ... |
| 17 | save_level1a | calibrationTool,logFile, calibratedSpectra, warningLevel0 | calibrationTool | Conditional | ... |

### Integration function
* read_meteo_data
* add_meteo_data
* read_level1a
* checking_channel_quality
* integrate_calibrated_spectra
* plot_integrated_spectra
* tropospheric_correction
* window_correction
* check_integrated
* save_level1b