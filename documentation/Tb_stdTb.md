# Brightness temperature computation

## Summary
Reminder note for the computation of the brightness temperature and its standard deviation.

## Tb
During a "standard" calibration cycle for GROSOM, we use the following procedure to compute the brightness temperature of the sky.

### Calibration

Here we assume that the calibration time is 10 minutes and that we got N measurement of the hot, cold and sky (in practice, N is different for each of the load).

1. Find the good hot and cold indices (see outlier detection)
2. Mean the hot and cold spectra for good indices as well as T_hot (recorded during the hot and cold measurement)
3. Find the good sky indices
4. Use calibration equation on ALL the good sky indices --> we get N brightness temperature spectra
5. Use the N spectra to compute the standard deviation of Tb at each channel during this calibration cycle. At this point, we also compute the mean std deviation for this calibration cycle.
6. Take the mean of the sky spectra corresponding to the good sky indices
7. Re-use the calibration formula with the mean sky spectrum for the calibration cycle. This is the one saved and integrated further !
8. In addition, we added the variable "noise_level" which tries to describe the amount of noise present in a given spectra. It is computed as the std deviation of the diff of the Tb spectrum. 

#### Calibration formula

Tb =  (T_hot - T_cold) * (U_sky - U_cold)/(U_hot - U_cold) + T_cold

### Integration

Consider we have M calibration cycle within each integration cycle.

During the integration sub-process, the following is done for the brightness temperature:

1. Find the good calibration cycle based on the calibration flags (see flags)

2. Integrated Tb is computed as the mean of the M calibrated spectra from each calibration cycle.

3. The std deviation of the integrated spectra is computed as:

sqrt(sum(variance(calibrated Tb))/M)

4. We then compute the mean of the spectrum of std Tb (mean_std_Tb in level1b)

5. In addition, we added the variable "noise_level" which tries to describe the amount of noise present in a given spectra. It is computed as the std deviation of the diff of the integrated Tb spectrum. 

The difference between the "noise_level" and the "mean_std_Tb" is that the formed express the temporal variation of the sky temperature during the integration time while the latter is a characterization of how noisy is the spectrum during this time.


### What happens if we compare a 10 min calibration time integrated on 1 hour and directly a 1 hour calibration time. 

We will get a difference !

There are multiple reasons behind this difference:
1. The mean hot and cold spectra as well as T_hot inserted in the calibration formula are different. As they appear at the denominator, it results in different Tb if the averaging of the hot (or cold) spectra is done after (mean of 6 10 min) or before (mean of 1h) insertion in the calibration formula.
2. Due to the nature of the outliers detection, the number of spectra averaged for hot, cold and sky might be different.
3. ...

The difference after retrieval might be interesting to do.

### Why we chose the 2 step procedure
Because of the variability of the atmosphere which can change on a time scale smaller than 1 hour. 