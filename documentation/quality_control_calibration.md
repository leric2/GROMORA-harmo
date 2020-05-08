# Quality control calibration

Summary of the quality control applied during the GROSOM calibration routine (level0 to level1b).

It can be divided into 2 main categories: the removing of "outlier" and the flagging of the data.

## Flags
I have tried to prefer the flagging of data instead of the removal of data whenever possible. During the calibration routine, I have set up the following flags to help understand why some data are missing (on a later stage) or showing some weird profiles after retrieval. 

There are different flags applied between the lvl0 to lvl1b routine.

#### Level0 flags:
When reading the raw data, we perform a few checks on the quality of it in the function *check_level0_generic()*:
a. size of log and raw data
b. check number of tipping curves done
c. check if the number of cycle is realistic

It writes a "warningLevel0" structure in the form of a string which is later saved in the lvl1a netCDF file as a global attributes. Note that we can deactivate this check by setting the "checkLevel0" parameters of the **calibrationTool** structure to False.

Additionaly, the "comment" variable present in the log file is also saved as a global attribues to the lvl1a file.

#### Calibration Level1a flags
These flags are defined into the "check_calibrated()" function after the effective calibration is done. The following is implemented:
a. "sufficientNumberOfIndices": checks that each calibration cycle averaged enough individual cycle of each type (hot, cold and antenna).
b. "systemTemperatureOK": check that the mean Tsys was in a certain interval during the calibration cycle and that its standard deviation was under a certain threshold during the cycle. 
c. "LN2SensorsOK": checks if there are some averaged individual cycles that had this flags in their log file. 
d. "LN2LevelOK": same as above but for the level of LN2 in the cold load.
e. "hotLoadOK": check that the hot load was in a certain interval during the calibration cycle and that the standard deviation of Thot was under a certain threshold during the cycle. 
f. "FFT_adc_overload_OK": checks that no overload of the FFT adc occured on any of the individual cycle averaged in this calibration cycle (flags in the log file already).

 
## Outlier detection
At some point of the routine, we are forced to remove some outliers data sothat they do not pollute the rest of the data for the day, the cycle, etc.. 

This is the case when performing the calibration with hot and cold load to avoid that spurious spectra takes too much power in the averaging on a cycle. 