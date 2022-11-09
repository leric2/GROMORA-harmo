# GROMORA level 0

## Summary
This file describes the raw file (level 0) from GROMOS and SOMORA. It is potentially applicable to other MW radiometers designed at the IAP as the fields entries are often quite similar.

The raw MW data actually contain two files:
* a binary file with the raw spectral spectrometers voltages or counts
* a text file (the so-called log file) which contains all additional meta data. Each line of this text file correspond to a single spectrum in the binary file.

Note that many parameters (especially the older analog sensors) were simply not used and some of them were measuring unknown quantities which are therefore not used in the calibration routine.

---


## Standard parameters

In order for the calibration routine to be as generic as possible, there are a set of standard log parameters in the log file. Most of them are present in the original log file from the instruments but most often with different naming conventions. Therefore, it was decided to harmonize the naming of the log variables so that it can be used later on in the routine in a generic way. In this section, we enumerate the different log variables needed to run the calibration and list the corresponding original (OG) variables names in GROMOS and SOMORA log files.

|variable | Description | GROMOS OG log name | SOMORA OG log name
|-------|------|:-----------|:-----------|
| file | filename | Same | Same
| comment | any comments for this log file | Same | Same
| x | matrix with all the log file numeric data | Same | Same
| header | header names | Same | Same
| Year | year (yyyy) | Same | Same
| Month | month (mm) | Same | Same
| Day | day (dd) | Same | Same
| Hour | hour (HH) | Same | Same
| Minute |  minute (MM)  | Same | Same
| Second |  second (ss)  | Same | Same
| time | time as datenum | None | None
| dateTime | time as datetime with timezone defined| None | None
| Position | the mirror position corresponding to the one defined in calibrationTool | Same | Same
| Elevation_Angle| elevation angle of the measurement | Same | Same
| T_Hot_Absorber | the temperature of the hot calibration target measured at the absorber (Kelvin) | T_Hot (in °C) | AI_0
| Tipping_Curve_active| 1 is this is a tipping curve cycle| Same | Position 
| T_Window | the window temperature in Kelvin | TExt3| AI_3
| Freq_Lock| 1 if the frequency is locked during the cycle | PLL_Lock & Ferranti_Lock | IF_LO1_Lock & IF_LO2_Lock
| FE_T_Sys| receiver noise temperature measured by labview for each calibration cycle | Same | Same
| LN2_Sensors_OK| 1 if LN2 sensor is ok | LN2_Relay | Same
| LN2_Level_OK| 1 if LN2 level is ok | LN2_above_High & LN2_above_Low | Same
| T_Room | room ambiant temperature | T_Ceiling | AI_1
| V_Gunn | Gunn voltage | Same| Same
| FFT_adc_range |range of the ADC from the FFT | Same| Same
| FFT_adc_overload |number of spectra with adc overload during this cycle | Same| Same
| FFT_T_FPGA | temperature of the spectrometer FPGA | Same| Same
| FFT_Mode | mode of the FFT ? | Same| Same
| FFT_Nr_of_acq | number of spectral acquisition for this cycle | Same| Same
| Spectr_left_wing_start | DEPRECATED | Same| Same
| Spectr_left_wing_width | DEPRECATED | Same| Same
| Spectr_line_center | DEPRECATED | Same| Same
| Spectr_line_width | DEPRECATED | Same| Same
| Spectr_T_Line_Amp | DEPRECATED | Same| Same
| Spectr_T_Peak | DEPRECATED | Same| Same
| Spectr_T_Wing | DEPRECATED | Same| Same
| Data_file_size | Number of the corresponding start byte in spectral binary file ? | Same| Same
| SW_version | raw software version (yyymmdd) | Same| Same
| t | time of day | Same| Same

---

---


## Instrument specific parameters

### GROMOS

For GROMOS, there are a few key changes made to the log file during the decade 2010-2020:
1. 10.03.2010: before this date, the positions of the mirror were coded as strings (e.g. "COLD")
2. 06.05.2010: the binary filename extension chaned from ***.dat*** to ***.bin***.
3. 13.05.2010: start of T_room measurement
4. 13.09.2011: introduction of new temperature sensors with much more measurements after this time (see below). 

For more information on the temporal changes of the log file, see the [harmonize_log_gromos](../scripts/calibration/harmonize_log_gromos.m) function in the calibration folder.

| variable | Description | Comments
|-------|------|:-----------|
| AI_0 | ? | |
| AI_1 | ? | |
| AI_2 | ? | |
| AI_3 | ? | |
| AI_4 | ? | |
| AI_5 | ? | |
| AI_6 | ? | |
| AI_7 | ceiling temperature | used before 13.09.11 |
| AI_8 | ? |  |
| AI_9 | ? | |
| AI_10 | ? |  |
| AI_11 | ? |  |
| AI_12 | ? |  |
| AI_13 | ? |  |
| AI_14 | ? |  |
| AI_15 | ? |  |
| DIO_0 | ? |  |
| DIO_1 | ? |  |
| DIO_2 | ? |  |
| Tipping_Angle_Nr | indice of the angle used when doing a tipping curve | not sure why this variable exist as the angle is taken from Elevation_Angle |
| Ferranti_Lock | 1 if Ferranti lock is ON|  |
| PLL_Lock | 1 if PLL lock is ON | |
| LN2_above_High | flag to indicate if LN2 level is too high | this variable and the next are very confusing and their value seemed to have change with time 
| LN2_above_Low |  flag to indicate if LN2 level is high enough |
| LN2_Relay | 1 if LN2 relay is not in use (0 is good) |
| TExt0 | T_Ceiling | since 13.09.11 |
| TExt1 | T_Floor |  since 13.09.11  |
| TExt2 | T_Aircon_Out |  since 13.09.11  |
| TExt3 | T_Window |  since 13.09.11  |
| TExt4 | T_Amp1 (first amplifier) |  since 13.09.11  |
| TExt5 | T_Amp2 (second amplifier) |   since 13.09.11 |
| TExt6 | T_Mirror_View |  since 13.09.11  |
| TExt7 | T_Reserved ? |  since 13.09.11  |
| IWV | estimation of the integrated water vapor | not sure this is computed |   

---

### SOMORA

SOMORA raw files are slightly more "stable" than the GROMOS ones with time.

For more information on the temporal changes of the log file, see the [harmonize_log_somora](../scripts/calibration/harmonize_log_somora.m) function in the calibration folder.


| variable | Description | Comments
|-------|------|:-----------|
| AI_0 | T_Hot_Absorber | |
| AI_1 | T_Room | |
| AI_2 | ? | |
| AI_3 | window temperature (T_Window) | |
| AI_4 | ? | |
| AI_5 | ? | |
| AI_6 | ? | |
| AI_7 | T_Out |outside temperature ? |
| AI_8 | ? |  |
| AI_9 | ? ||
| AI_10 | ? | |
| AI_11 | ? | |
| AI_12 | ? | |
| AI_13 | ? | |
| AI_14 | ? | |
| AI_15 | ? | |
| DIO_0 | ? | |
| DIO_1 | ? | |
| DIO_2 | ? | |
| T_Hot | hot load temperature | This one is located at the blower and was off for almost 2 years (2015-2017) |
|IF_LO1_Lock  | lock of the first intermediate frequency | |
|IF_LO2_Lock  | lock of the first intermediate frequency | |