# GROSOM Meeting

23.07.2020

## Calibration
* Operationnal for SOMORA after 2013 and for GROMOS after 2018 (meteo data missing) (see demo folder)

### Quality control 
Open to discussion, done as follow right now:

#### Level0-1a
1. During effective calibration, we check individual spectra for:

    1. Elevation angle check for hot, cold and sky position. Absolute for hot and cold (might change for some period) and with a threshold (5deg) for sky measurement.
    
    2. FFTC_adc_ovedrload in log file (was done globally and might be too much for some days where overload was constant, but then should we keep these days ?)

    3. For hot and cold counts, we check that every cycle has at least 5% of its channels counts within the daily median counts +- 3 std deviation. (could be made on the 10 minutes averages but works quite good like that) We do the same for sky measurementt with +- 6 std deviation around the 10 minutes median.

All individual cycle removal is recorded and plotted (if not too many) but not saved. 

At that point, we also use the indiviual (cleaned) cycle to compute the standard deviation of TSys and Tb, using the global mean of all channel for Tsys and keeping a value for each channel in Tb.

2. After the calibration, we have a set of checks that are runned on the calibration cycle and are setting the following flags:

  		calibration_flags:errorCode_1 = "sufficientNumberOfIndices" ;
  		calibration_flags:errorCode_2 = "systemTemperatureOK" ;
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

#### Level1a-1b
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
3. raw_file_warning: 3 tyoes of warning can be found here (as string):
    1. consistency:unconsistentBinarySize
    2. consistency:numberCycles
    3. extra_timestamp

### Work to do next on calibration side:
1. Decide on quality control and adapt the routine until the beginning of FFTS
2. Final work of code documentation when the routine is complete
3. Identification of level 1 suspect data and link it with the labview log.
9. Check and assess the calibration quality ?

## Retrieval

Basic structure and workflow seem to work (at least for 2 different days and both instrument).
For now it works with the following parameter:
* ARTS development version n. 2.3.1285
* tropospheric corrected spectra as measurement (Ingold, v1)
* Fixed noise value 0.1K, no binning but bad channels identification and noise value raised for "bad channels"
* PTZ profiles from ECMWF oper dataset and CIRA86. For now, using only the closest lat-lon and 6h time resolution.
* apriori ozone profiles from old GROMOS and SOMORA retrievals respectively
* apriori cov for ozone set to constant 2ppm
* no baseline, no freq shift, ...
* Spectroscopy from Perrin (default coming with the new ARTS)
* Species defined as in GROMOS old retrieval ["O3","H2O-PWR98", "O2-PWR93","N2-SelfContStandardType"], to discuss ? it will also need to be changed for including water vapor.
* FM: simulation grid from 800 to 114e3m
* retrieval grid from 1000 to 95e3m

### Work to do next on retrieval side:
1. Work on understanding ARTS and its possibilities a bit more, especially to succed in point 2...
2. Include water vapor retrievals and use the non-tropospheric corrected spectra
3. Use pyarts instead of Tython in pyretrievals
4. Make a proper sensors definition and check how to include the DSB correction in ARTS
5. Decide which apriori we want to use:
    * O3: GEOMS ?
    * H2O: ECMWF + Fascod
6. Investigate properly the covariance matrices for measurement and apriori profiles
7. Spectroscopy to check, use HITRAN ?
8. Think about the outputs for level2
9. Check and assess the retrieval quality ?

## Atmospheric data
As atmospheric data, I am now using a "simplified" algorithm to extract the atmospheric state from the ECMWF oper dataset (6h resolution) at Payerne and Bern (closest grid point). This can be improve but is maybe not a priority at the moment.

General workflow is:
* Download ECMWF global file data (should we go to ERA5 then ?)
* Pre-processing of ECMWF data with pyretrievals for each locations

Then, an idea would be to store the apriori dataset separated from the level2 files and read it when needed (-> in the form of the ECWF store). Might be better than "creating" the ptz profile at each retrieval ? 

## Timeline
1. Holiday next week and early September (1 week)
2. TBA, depending on the priorities.

