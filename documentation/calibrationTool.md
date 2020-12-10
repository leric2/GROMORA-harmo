# calibrationTool structure

## Summary
This is the structure containing all the information to launch a calibration of GROSOM intrument. It contains not only all parameters required for calibrating a given instrument or a given day
but also the functions that will be used to performed the calibration.

All parameters must contain: type, units, where it is defined

## Parameters


| # | name | type | unit | Description |
|---|------|------|------|:-----------:|
| 1 | dateStr | str | unit | 'YYYY_MM_HH' |
| 2 | Year | str | double | YYYY |
| 3 | calibratedSpectraSpectralPlot | boolean | double | YYYY |

This is a table listing all the parameter stored into calibrationTool.

### Variable list
* dateStr (str):        'YYYY_MM_HH'
* Year (double):         YYYY
* Month
* Day
* dateTime
* timeNumber
* instrumentName
* lightSpeed
* h
* kb
* zeroDegInKelvin
* backgroundMWTb
* meanDatetimeUnit
* referenceTime
* calendar
* hotSpectraNumberOfStdDev
* coldSpectraNumberOfStdDev
* labviewLog
* calType
* rawSpectraPlot
* calibratedSpectraPlot
* calibratedSpectraSpectralPlot
* calibratedSpectraStdTbPlot
* integratedSpectraPlot
* calibrationTime
* integrationTime
* minNumberOfAvgSpectra
* filterByTransmittance
* filterByFlags
* TCold
* dataLocation
* PI_NAME
* PI_AFFILIATION
* PI_ADDRESS
* PI_EMAIL
* dataSource
* lon
* lat
* altitude
* azimuthAngle
* timeZone
* observationFreq
* numberOfSpectrometer
* spectrometer
* samplingRateFFTS
* numberOfChannels
* IQProcessing
* fLO1
* fLO2
* LOFreqTot
* instrumentBandwidth
* badChannels
* numberOfAquisitionSpectraHot
* numberOfAquisitionSpectraAntenna
* numberOfAquisitionSpectraCold
* binaryDataExtension
* logFileDataExtension
* bytesPerValue
* binaryType
* rawFileFolder
* extraFileFolder
* level1Folder
* filename
* file
* checkLevel0
* delimiter_logfile
* THotUnit
* harmonize_log
* positionIndAsName
* indiceCold
* indiceAntenna
* indiceHot
* elevationAngleAntenna
* elevationAngleCold
* elevationAngleHot
* elevationAngleTolerance
* elevationAngleHotTol
* elevationAngleColdTol
* cycleDurationCold
* cycleDurationSky
* cycleDurationHot
* flipped_spectra
* flip_spectra
* numberOfTippingCurveExpected
* toleranceTippingCurves
* goodFlagLN2Above
* goodFlagLN2Below
* numberOfCyclesExpected
* toleranceNumberCycles
* tippingSize
* TSysCenterTh
* TSysThresh
* stdTSysThresh
* THotTh
* THotAbsThresh
* hotTemperatureStdThreshold
* stdAntAngleThresh
* minNumberOfIndicePerCycle
* threshNumRawSpectraHot
* threshNumRawSpectraCold
* threshNumRawSpectraAnt
* maxProportionOfIndLN2LevelOutlier
* maxProportionOfIndLN2SensorOutlier
* frequencyBandAroundCenterTSys
* maxStdDevTbCal
* maxStdDevTbInt
* filterTypeChannelQualityCal
* filterTypeChannelQualityInt
* filter1
* filter2
* meteoFolder
* read_meteo_data
* add_meteo_data
* doTippingCurve
* tWindow
* troposphericCorrection

* filenameLevel1a
* successfulCalibration
* flagVectorLength
* logFile
* filenameLevel1b
* successfulIntegration

## Functions

* read_level0
* check_level0
* reformat_spectra
* plot_raw_spectra
* calibrate
* plot_calibrated_spectra
* check_calibrated
* save_level1a
* read_level1a
* checking_channel_quality
* integrate_calibrated_spectra
* plot_integrated_spectra
* tropospheric_correction
* window_correction
* check_integrated
* save_level1b