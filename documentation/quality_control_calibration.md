# Quality control calibration

## Summary

Technical documentation concerning the quality control applied during the GROSOM
calibration routine (level 0 to level 1b). The general principles of flags and
outlier detection is presented in the GROSOM user guide (UG) that can be found in the [User Guide](https://git.iap.unibe.ch/IAP_MCH/UserGuideGROSOM) repository. 

In general, the quality control of the calibration can be divided into 2 main
categories: the removing of "outlier" and the flagging of the potential spurious
data. The flagging of data is always preferred however, it is sometimes
necessary to remove some strong outliers that would otherwise "pollute" some good data (for instance when averaging).

The flags are indicative of the supposed quality of some data and can be used mainly in 2 ways:
1. To filter the "good" or "bad" time periods (e.g. during integration) to keep only good quality data in an operational retrieval for instance.
2. To understand the origin of some weird results arising in some time periods. 

There are different flags and outlier detection steps applied between the level 0 to level 1b routine.

## Level0 flags:

The first flags with a completely indicative nature, are the level 0 warnings.

When reading the raw data, we perform a few checks on their quality in the function *check_level0*:
a. size of log and raw data
b. check number of tipping curves done
c. check if the number of cycle is realistic

It writes a "warningLevel0" structure in the form of a string which is later saved in the lvl1a netCDF file as a global attributes. Note that we can deactivate this check by setting the "checkLevel0" parameters of the **calibrationTool** structure to False.

Additionaly, the "comment" variable present in the log file is also saved as a global attribues to the lvl1a file.

Saves extra timestamps in extra raw and log files.

## Calibration sub-routine

### Outlier detection during the calibration

At some point of the routine, we are forced to remove some outliers data sothat
they do not pollute the rest of the data for the day, the cycle, etc.. This is
the case when performing the calibration with hot and cold load to avoid that
spurious spectra takes too much power in the averaging on a calibration cycle. 

During the calibration process (in the *calibrate* function), we therefore apply some outliers detection to
identify and remove spurious individual hot, cold and sky spectra. More
specifically, we introduce the following outliers detection techniques or
variables to identify spurious spectra:

During effective calibration, we check individual spectra for:

1. Elevation angle check for hot, cold and sky position with an expected value and a certain threshold.
    
2. State of the analog-to-digital converter when a the spectrum was recorded (*FFT_adc_overload* in the log file)

3. Check that enough individual channels are comprised within a certain interval from the daily median (for hot and
        cold spectra). For sky spectra, as the atmosphere has a natural variation, we use the median of the calibration
        cycle instead of a daily median.

All the parameters controllings the outlier detecion are defined within the calibrationTool structure and can be adapted easily if needed. There is also a possibility to change the outlier detection techniques applied during a certain calibration by using the *outlierDectectionType* option in *calibrationTool*. It takes the following values:
* standard: all 3 detection techniques are applied
* noFFT: only 1. and 2. above, no check of the analog-to-digital converter applied (we realise that some time periods actually have constantly the *FFT_adc_overload* flag on).
* none: no outlier detection during the calibration

Note that individual cycle removal is recorded and plotted (except if this there are too many) but not saved. These are the red, blue or green crosses on the level 1a plots (see UG) to represents respectively hot, cold or sky outliers detected during a calibration cycle.

At that point, we also use the indiviual (cleaned) cycle to compute the standard deviation of the noise receiver temperature (stdTNoise) and of the brightness temperature (stdTb), using the channel averaged value of TNoise.

### Some possible Improvements

Add threshold on the number of ADC overloads per cycle ?

### Level 1a flags

In addition to the outliers removal during the calibration, there are some further checks done on the calibrated data.
These checks are done within the *check_calibrated* function and result in the creation of a set of
flags that are saved in the level 1a. Contrary to the outliers detection, the flags are indicative of the data quality
and do not lead to any data removal before level1a. Depending on the final use of the data, the user is then free to
take these flags into account of not. Also the flags are determined for each calibration cycle (and not on individual
spectrum) and are concerned with the followings parameters:

After the calibration, we have a set of checks that are runned on the calibration cycle and are setting the following flags:

  		calibration_flags:errorCode_1 = "sufficientNumberOfIndices" ;
  		calibration_flags:errorCode_2 = "noiseTemperatureOK" ;
  		calibration_flags:errorCode_3 = "LN2SensorsOK" ;
  		calibration_flags:errorCode_4 = "LN2LevelOK" ;
  		calibration_flags:errorCode_5 = "hotLoadOK" ;
  		calibration_flags:errorCode_6 = "PointingAngleOK" ;

1. check that a minimum of indices has been averaged for each position within the calibration.
2. For checking Tsys, it is open for discussion but I am doing the following:
    * Computing Tsys at the center channels (+- 200MHz) while removing some potential bad channels and comparingthis value to a threshold.
    * Using the std dev of Tsys computer during the calibration against a threshold.
3. Check that 30% of the LN2 Sensors flags in the log are OK
4. Check that 30% of the LN2 Level flags in the log are OK
5. Check that the std dev of the hot load temperature is below a certain threshold during the calibration cycle
6. Check that the std dev of the pointing angle of sky measurement is below a certain threshold during thecalibration cycle

## Integration sub-routine
During the second step of the routine, the following quality control is performed:
1. We check the channel quality on calibration cycle before doing a window and a generic tropospheric correction on each calibration cycle (not taking into accound the bad channels in the wings averaging for instance) --> this is then used for the integration.

2. During the integration, we average only the calibrated 10 minutes spectra that have:
    1. No calibration flags
    2. A tropospheric transmittance > 0.2 
This point is opened to discussion...

3. After the integration, we perform one more time a channel quality check, then perform a window and a troposheric correction for integrated (1h) cycle. We then output a 2 flags vector based on:

  		calibration_flags:errorCode_1 = "sufficientNumberOfAvgSpectra" ;
  		calibration_flags:errorCode_2 = "rain_Accumulation_OK" ;

    1. We check that the number of averaged spectra was at least 3 (30 minutes) which means it also contains a minimum of each individual cycles. 
    2. We check that the amount of rain was not too hight during this integration cycle. UNIT ??


        calibrationTool.filterByTransmittance = true;
        calibrationTool.transmittanceThreshold = 0.2;
        calibrationTool.filterByFlags = true;
        
        filtering_of_calibrated_spectra

## How to add flags or variables ?

#### Additional warnings as attributes of the level1a-b file
1. labview log warning: check if an entry is present in the labview log file for the day
2. raw_file_comment: if a comment is found in the raw file
3. raw_file_warning: 3 types of warning can be found here (as string):
    1. consistency:unconsistentBinarySize
    2. consistency:numberCycles
    3. extra_timestamp

