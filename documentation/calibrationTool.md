# calibrationTool structure

## Summary
This is the structure containing all the information to launch a calibration of GROSOM intrument. It contains not only all parameters required for calibrating a given instrument or a given day
but also the functions that will be used to performed the calibration.

All parameters must contain: type, units, where it is defined

## Table of Contents
1. [Parameters](#Parameters)
2. [Functions](#Functions)

## Parameters

### Time variable (11)

|variable | type  | Description |
|------|------|:-----------:|
| dateStr | str  | 'YYYY_MM_HH' |
| Year | double  | YYYY |
| Month | double  | MM |
| Day | double | DD |
| [dateTime](dateTime) | datetime |  |  
| timeNumber | datenum |  |  
| meanDatetimeUnit |  |  |  
| referenceTime |  |  |  
| calendar |  |  |  
| timeZone |  |  |  
| calendar |  |  |  

#### dateTime

A matlab structure datetime. Defined with timeZone

---

### Metadata (7)
|variable | type | Description |
|-------|------|:-----------:|
| instrumentName | str |  |  
| dataLocation |  |  | 
| PI_NAME |  |  | 
| PI_AFFILIATION |  |  | 
| PI_ADDRESS |  |  | 
| PI_EMAIL |  |  | 
| dataSource |  |  | 


---


### Geolocation data (4)
|variable | type | Description |
|------|------|:-----------:|
| lon |  |  | 
| lat |  |  | 
| altitude |  |  | 
| azimuthAngle |  |  | 


---


### Physical constant / parameters (6)

| variable | type | unit | Description |
|------|------|------|:-----------:|
| lightSpeed | double | m/s |  |
| h | double |  |  |
| kb | double |  |  |
| zeroDegInKelvin | double |  |  |
| backgroundMWTb | double |  | MW background radiation |

---

### Spectrometer variables (15-16)
| variable | type | Description |
|------|------|:-----------:|
| observationFreq | double |  |  
| numberOfSpectrometer | double |  |  
 | spectrometer |  |  |  
 | samplingRateFFTS | double |  |  
 | numberOfChannels | double |  |  
 | IQProcessing |  |  |  
 | fLO1, fLO2, ... | double |  |  
 | LOFreqTot | double |  |  
 | instrumentBandwidth | double | bandwidth \[Hz\] |  
 | badChannels |  |  |  
 | numberOfAquisitionSpectraHot |  |  |  
 | numberOfAquisitionSpectraAntenna |  |  |  
 | numberOfAquisitionSpectraCold |  |  |  
 | flipped_spectra | boolean  |  |  

---

### Raw files check (6)

| variable | type | Description |
|------|------|:-----------:|
| checkLevel0 | boolean |  | 
| tippingSize | double |  |  
| numberOfTippingCurveExpected | double |  | 
| toleranceTippingCurves | double |  | 
| numberOfCyclesExpected | double |  | 
| toleranceNumberCycles | double |  | 

---

### Files and folder variables (14)
| variable | type  | Description |
|------|------|:-----------:|
| binaryDataExtension | str |  | 
| logFileDataExtension |  |  |  
| bytesPerValue |  |  |  
| binaryType |  |  |  
| rawFileFolder |  |  |  
| extraFileFolder |  |  |  
| level1Folder |  |  |  
| meteoFolder | str |  | 
| filename |  |  |  
| [file](#112) |  |  |  
| delimiter_logfile |  |  |  
| labviewLog | boolean |  |  
| filenameLevel1a |  |  |  
| filenameLevel1b |  |  |  

#### 112


---

### Log variables (14)

| variable | type  | Description |
|------|------|:-----------:|
| THotUnit | str |  |  
| positionIndAsName | boolean |  |  
| indiceCold | double |  |  
| indiceAntenna | double |  |  
| indiceHot | double |  |  
| elevationAngleCold | double |  |  
| elevationAngleAntenna | double |  |  
| elevationAngleHot | double |  |  
| elevationAngleColdTol | double |  | 
| elevationAngleTolerance | double |  | 
| elevationAngleHotTol | double |  | 
| cycleDurationCold | double |  | 
| cycleDurationSky | double |  | 
| cycleDurationHot | double |  | 

---

### Calibration variables + outlier detection ? (22)
| variable | type  | Description |
|------|------|:-----------:|
| calType | str |  | 
| calibrationTime |  | s |  
| TSysCenterTh |  |  |  
| TSysThresh |  |  |  
| stdTSysThresh |  |  |  
| frequencyBandAroundCenterTSys |  |  |  
| THotTh |  |  |  
| THotAbsThresh |  |  |  
| hotTemperatureStdThreshold |  |  |  
| stdAntAngleThresh |  |  |  
| minNumberOfIndicePerCycle |  |  |  
| hotSpectraNumberOfStdDev |  |  |  
| coldSpectraNumberOfStdDev |  |  |  
| threshNumRawSpectraHot |  |  |  
| threshNumRawSpectraCold |  |  |  
| threshNumRawSpectraAnt |  |  | 
| maxProportionOfIndLN2LevelOutlier |  |  | 
| maxProportionOfIndLN2SensorOutlier |  |  | 
| threshNumRawSpectraAnt |  |  |
| maxStdDevTbCal |  |  |  
| goodFlagLN2Above |  |  |  
| goodFlagLN2Below |  |  |  

---

### Meteo variables (3)
| variable | type  | Description |
|------|------|:-----------:|

* doTippingCurve
* tWindow
* troposphericCorrection

---

### Plot variables  (4)

| variable | type  | Description |
|---|------|------|:-----------:|
| rawSpectraPlot | boolean |  | 
| calibratedSpectraPlot | boolean |  | 
| calibratedSpectraSpectralPlot | boolean |  | 
| integratedSpectraPlot | boolean |  | 

---

### Integration variables (9)
| variable | type  | Description |
|---|------|------|:-----------:|
| integrationTime |  | s |  
| minNumberOfAvgSpectra | double |  |  
| filterByTransmittance | boolean |  |  
| filterByFlags | boolean |  |  
| filterTypeChannelQualityCal | int |  |  
| filterTypeChannelQualityInt | int |  |  
| filter1 | struct |  |  
| filter2 | struct |  |  
| maxStdDevTbInt | boolean |  | 

---

### Variable list
* successfulCalibration
* flagVectorLength
* logFile
* 
* successfulIntegration

---
---

## Functions

Functions stored in calibrationTool only.

### Calibration function

All functions here are used in [run_calibration](run_calibration.md) and are detailed there.

| name | inputs | outputs | type | Description |
|------|------|------|------|:-----------:|
| read_level0 | calibrationTool, readingRawFile | logFile, rawSpectra | Required | import raw files |
| harmonize_log | calibrationTool, logFile | logFile | Required |
| reformat_spectra | rawSpectra, logFile, calibrationTool| rawSpectra | Conditional |
| check_level0 | logFile, rawSpectra, calibrationTool | warningLevel0 | Optional |
| flip_spectra | rawSpectra | rawSpectra | Conditional |
| plot_raw_spectra | rawSpectra, plot parameters | - | Optional |
| read_meteo_data | calibrationTool | logFile.meteo | Optional ?|
| run_tipping_curve | rawSpectra, logFile, calibrationTool | logFile.TC | Conditional |
| calibrate | rawSpectra, logFile, calibrationTool, calibrationTool.calType | drift, calibratedSpectra | Required |
| check_calibrated | logFile, calibrationTool,calibratedSpectra | calibratedSpectra | Conditional |
| plot_calibrated_spectra | calibrationTool, drift, logFile.meteo, calibratedSpectra, N | - | Conditional |
| save_level1a | calibrationTool,logFile, calibratedSpectra, warningLevel0 | calibrationTool | Conditional |

---

### Integration function

All functions here are used in [run_integration](run_integration.md) and are detailed there.

| name | inputs | outputs | type | Description |
|------|------|------|------|:-----------:|
| read_level1a |calibrationTool |calibratedSpectra, meteoData, calibrationTool | Required |
| add_meteo_data | calibrationTool, meteoData, calibratedSpectra | calibratedSpectra | Required |
| check_channel_quality | calibratedSpectra,calibrationTool,filterType | calibratedSpectra | Required |
| tropospheric_correction | calibratedSpectra, calibrationTool | calibratedSpectra | Required |
| integrate_calibrated_spectra | calibrationTool,calibratedSpectra | integratedSpectra | Required |
| check_integrated | calibrationTool, integratedSpectra | integratedSpectra | Required |
| plot_integrated_spectra | calibrationTool, integratedSpectra | - | Required |
| save_level1b | calibrationTool, level1 | calibrationTool | Required |