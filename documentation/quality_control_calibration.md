# Quality control calibration

## Summary

Technical documentation concerning the quality control applied during the GROMORA
calibration routine (level 0 to level 1b). The general principles of flags and
outlier detection is presented in the GROMORA user guide (UG) that can be found in the [User Guide](https://git.iap.unibe.ch/IAP_MCH/UserGuideGROMORA.git) repository. 

In general, the quality control of the calibration can be divided into 2 main
categories: the removing of "outlier" and the flagging of the potential spurious
data. The flagging of data is always preferred however, it is sometimes
necessary to remove some strong outliers that would otherwise "pollute" some good data (for instance when averaging).

The flags are indicative of the supposed quality of some data and can be used mainly in 2 ways:
1. To filter the "good" or "bad" time periods (e.g. during integration) to keep only good quality data in an operational retrieval for instance.
2. To understand the origin of some weird results arising in some time periods. 

There are different flags and outlier detection steps applied between the level 0 to level 1b routine.

## Labview log warning

A very first basic check is made by reading the electronic logbook of the raw
labview software. In case an entry exist for the day on which the calibration is
performed, this will be displayed in the global attribute
(*labview_logfile_warning*) in the level 1a and 1b which can take the following values:
* "check labview log !": means that an entry is recorded in the labview log for
  this day
* "clean": means that no entry was found for this day.

## Level0 flags:

The first flags with a completely indicative nature, are the level 0 warnings.

When reading the raw data, we perform a few checks on their quality in the function *check_level0*:
1. size of log and raw data
1. check number of tipping curves done
1. check if the number of cycle is realistic
1. extra_timestamp: indicates if some additional time stamps corresponding to
   another day were found in the log file. If this is the case, these have been
   saved in extra raw and log files that can be read with the corresponding day.

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

1. Elevation angle check for hot, cold and sky position with an expected value
   and a certain threshold.

        calibrationTool.elevationAngleCold
        calibrationTool.elevationAngleAntenna
        calibrationTool.elevationAngleHot
        calibrationTool.elevationAngleColdTol
        calibrationTool.elevationAngleTolerance
        calibrationTool.elevationAngleHotTol
    
2. State of the analog-to-digital converter when a the spectrum was recorded
   (*FFT_adc_overload* in the log file). We check that the number of ADC
   overloads per cycle does not exceed a certain threshold.

        calibrationTool.adcOverloadThresh 


3. Check that enough individual channels are comprised within a certain interval from the daily median (for hot and
        cold spectra). For sky spectra, as the atmosphere has a natural variation, we use the median of the calibration
        cycle instead of a daily median.

        calibrationTool.threshNumRawSpectraHot
        calibrationTool.hotSpectraNumberOfStdDev
        calibrationTool.threshNumRawSpectraCold
        calibrationTool.coldSpectraNumberOfStdDev
        calibrationTool.threshNumRawSpectraAnt
        calibrationTool.skySpectraNumberOfStdDev

All the parameters controlling the outlier detecion are defined within the
*calibrationTool* structure (shown above) and can be adapted easily if needed. There is
also a possibility to change the outlier detection techniques applied during a
certain calibration by using the *outlierDectectionType* option in
*calibrationTool*. It takes the following values:
* standard: all 3 detection techniques are applied
* noFFT: only 1. and 2. above, no check of the analog-to-digital converter applied (we realise that some time periods actually have constantly the *FFT_adc_overload* flag on).
* none: no outlier detection and removal applied during the calibration

Note that individual cycle removal is recorded and plotted (except if this there are too many) but not saved. These are the red, blue or green crosses on the level 1a plots (see UG) to represents respectively hot, cold or sky outliers detected during a calibration cycle.

At that point, we also use the indiviual (cleaned) cycle to compute the standard deviation of the noise receiver temperature (stdTNoise) and of the brightness temperature (stdTb), using the channel averaged value of TNoise.

### Some possible Improvements
 

### Level 1a flags

In addition to the outliers removal during the calibration, there are some
further checks done on the calibrated data. These checks are done within the
*check_calibrated* function and result in the creation of a set of flags that
are saved in the *flags* group of level 1a (see [level 1](level1.md)). Contrary
to the outliers detection, the flags are indicative of the data quality and do
not lead to any data removal before saving level 1a. Depending on the final use of the
data, the user is then free to take these flags into account of not. Also the
flags are determined for each calibration cycle (and not on individual
spectrum). Note that all the below listed flags can be changed or setup in the
*calibrationTool* structure as shown for each flag.

1. Check that a minimum of indices has been averaged for each position within
   the calibration.

  		calibrationTool.minNumberOfIndicePerCycle

2. The noise receiver temperature (*TNoise*): 

Compute *TNoise* at the center channels (standard value is +- 200MHz around line
   center) while removing some potential bad channels and compare it with a
   standard value with a threshold. 
   
  		calibrationTool.frequencyBandAroundCenterTNoise
  		calibrationTool.TNoiseCenterTh 
        calibrationTool.TNoiseThresh

Compare the std dev of *TNoise* computed during the calibration (from the drift
   structure) against a threshold.

        calibrationTool.stdTNoiseThresh

3. Check that X% (standard is 30%) of the LN2 Sensors flags in the log are OK
   for this calibration cycle.

        calibrationTool.maxProportionOfIndLN2SensorOutlier

4. Check that X% (standard is 30%) of the LN2 Level flags in the log are OK for
   this calibration cycle.

        calibrationTool.maxProportionOfIndLN2LevelOutlier

5. Check that the standard deviation of the hot load temperature is below a
   certain threshold during the calibration cycle

        calibrationTool.hotTemperatureStdThreshold

6. Check that the standard deviation of the pointing angle of sky measurement is below a
   certain threshold during the calibration cycle

        calibrationTool.stdAntAngleThresh

The flags are then saved in the form of a boolean vector in the level 1a netCDF
files and the correspondance of each vector element to each flag is added as
attributes:

  		calibration_flags:errorCode_1 = "sufficientNumberOfIndices" ;
  		calibration_flags:errorCode_2 = "noiseTemperatureOK" ;
  		calibration_flags:errorCode_3 = "LN2SensorsOK" ;
  		calibration_flags:errorCode_4 = "LN2LevelOK" ;
  		calibration_flags:errorCode_5 = "hotLoadOK" ;
  		calibration_flags:errorCode_6 = "PointingAngleOK" ;

Note that a full vector of 1 is indicative of a good quality data, while the
presence of 0 indicates a problem within the calibration cycle.

During the calibration sub-routine, any flagged cycle is displayed with the
corresponding description of the flag(s). In addition, the flagged calibrated cycles are
identified as magenta crosses on the level 1a plots (see UG).

## Integration sub-routine

During the second step of the routine, different kind of quality control are
made, first on the calibrated spectra read from the level 1a and then on the
integrated spectra. 

### Channel quality

Some spurious channels are often present in the spectrometer data and it was
decided to keep all channels as long as possible during the processing to avoid
loosing data. This implies however, that the channel quality must be verified
before doing any channel averaging if we want to avoid spurious channels to
impact too much the results. 

With that in mind, we check the channel quality of the calibrated spectra before
doing the window and tropospheric correction which are then used for the
selection of data for integration. 

### Integration of the calibrated spectra

When integrating together the calibrated spectra, there are different options
that can be chosen depending on the type of analysis to perform later.
Therefore, we introduced some simple filtering techniques to choose which
calibrated spectra to keep within an integration cycle. For all filters, it is
possible to activate and/or change the thresholds (of the transmittance value
for instance) in the *calibrationTool* structure. For now, we have 2 filters
that are concerned with:

1. The flags of the calibrated spectra.

        calibrationTool.filterByFlags = true;

2. The value of the tropospheric transmittance from each calibrated spectra. 

        calibrationTool.filterByTransmittance = true;
        calibrationTool.transmittanceThreshold = 0.2;

The type of filtering used during the integration is saved in the
level 1b as the *filtering_of_calibrated_spectra*. global attributes.

During the calibration sub-routine, any flagged cycle is displayed with the
corresponding description of the flag(s). In addition, the flagged integrated cycles are
identified as magenta stars on the level 1b plots (see UG).

### Level 1b flags
After the integration, we perform one more time a channel quality check, then
perform a window and a troposheric correction for the integrated (1h) spectra.
Similar to the calibration sub-routine, we introduce a set of flags based on:
1. The number of averaged calibrated spectra within this integrated spectra
   (which means it also contains a minimum of each individual cycles). 
2. The value of the transmittance for the integrated spectra is higher than a
   certain threshold (note that the threshold can be different than the ones
   used for the filtering of the calibrated spectra)

        calibrationTool.troposphericTransmittanceFlag

It results in the creation of a boolean vector which is saved in the level 1b,
also in the *flags* group with the corresponding attributes:

  		calibration_flags:errorCode_1 = "sufficientNumberOfAvgSpectra" ;
  		calibration_flags:errorCode_2 = "tropospheric_transmittance_OK" ;

### Some possible Improvements

Consider only certain flags during the integration to select the "good"
calibrated spectra (in case *filterByFlags* is true)

## How to add flags ?

When the time of reprocessing will come or in the case we want to adapt the
routine to a new instrument, this will most likely be needed to adapt the flags
for level 1a or 1b. This needs to deep a bit into the codes but should be quite
easy to do.

1. Flags are defined in the *check_calibrated* (*check_integrated*) functions
   for the level 1a (level 1b). 

2. New flags should be defined to be 1 when OK and 0 when not. 

3. New flags can be added at the end of the *errorVector* while the description
   of the new flags should be added at the end of *errorVectorDescription*.

4. The saving of the flags and their description should be done automatically
   but it might be worth checking ;)



