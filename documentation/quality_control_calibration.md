# Quality control calibration

Summary of the quality control applied during the GROSOM calibration routine (level0 to level1b).

It can be divided into 2 main categories: the removing of "outlier" and the flagging of the data.


I have tried to prefer the flagging of data instead of the removal of data whenever possible. During the calibration routine, I have set up the following flags to help understand why some data are missing (on a later stage) or showing some weird profiles after retrieval. 

There are different flags applied between the lvl0 to lvl1b routine.

## Level0 flags:
When reading the raw data, we perform a few checks on the quality of it in the function *check_level0_generic()*:
a. size of log and raw data
b. check number of tipping curves done
c. check if the number of cycle is realistic

It writes a "warningLevel0" structure in the form of a string which is later saved in the lvl1a netCDF file as a global attributes. Note that we can deactivate this check by setting the "checkLevel0" parameters of the **calibrationTool** structure to False.

Additionaly, the "comment" variable present in the log file is also saved as a global attribues to the lvl1a file.

Saves extra timestamps in extra raw and log files.

## Outlier detection during the calibration
At some point of the routine, we are forced to remove some outliers data sothat they do not pollute the rest of the data for the day, the cycle, etc.. 

This is the case when performing the calibration with hot and cold load to avoid that spurious spectra takes too much power in the averaging on a cycle. 

1. During effective calibration, we check individual spectra for:

    1. Elevation angle check for hot, cold and sky position. Absolute for hot and cold (might change for some period) and with a threshold (5deg) for sky measurement.
    
    2. FFTC_adc_ovedrload in log file (was done globally and might be too much for some days where overload was constant, but then should we keep these days ?)

    3. For hot and cold counts, we check that every cycle has at least 5% of its channels counts within the daily median counts +- 3 std deviation. (could be made on the 10 minutes averages but works quite good like that) We do the same for sky measurementt with +- 6 std deviation around the 10 minutes median.

All individual cycle removal is recorded and plotted (if not too many) but not saved. 

At that point, we also use the indiviual (cleaned) cycle to compute the standard deviation of TSys and Tb, using the global mean of all channel for Tsys and keeping a value for each channel in Tb.

## Calibration Level1a flags

2. After the calibration, we have a set of checks that are runned on the calibration cycle and are setting the following flags:

  		calibration_flags:errorCode_1 = "sufficientNumberOfIndices" ;
  		calibration_flags:errorCode_2 = "noiseTemperatureOK" ;
  		calibration_flags:errorCode_3 = "LN2SensorsOK" ;
  		calibration_flags:errorCode_4 = "LN2LevelOK" ;
  		calibration_flags:errorCode_5 = "hotLoadOK" ;
  		calibration_flags:errorCode_6 = "PointingAngleOK" ;

    1. check that a minimum of indices has been averaged for each position within the calibration.
    2. For checking Tsys, it is open for discussion but I am doing the following:
        * Computing Tsys at the center channels (+- 200MHz) while removing some potential bad channels and comparing this value to a threshold.
        * Using the std dev of Tsys computer during the calibration against a threshold.
    3. Check that 30% of the LN2 Sensors flags in the log are OK
    4. Check that 30% of the LN2 Level flags in the log are OK
    5. Check that the std dev of the hot load temperature is below a certain threshold during the calibration cycle
    6. Check that the std dev of the pointing angle of sky measurement is below a certain threshold during the calibration cycle

## Level1a-1b
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


#### Additional warnings as attributes of the level1a-b file
1. labview log warning: check if an entry is present in the labview log file for the day
2. raw_file_comment: if a comment is found in the raw file
3. raw_file_warning: 3 types of warning can be found here (as string):
    1. consistency:unconsistentBinarySize
    2. consistency:numberCycles
    3. extra_timestamp
