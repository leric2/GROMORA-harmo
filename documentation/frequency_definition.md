# Frequency Vector Definition for FFTS
## Summary
As have been discussed multiple time, the definition of the frequency vector corresponding to the FFTS channel is not a straightforward task and this was done differently for GROMOS and SOMORA in the past. We will try to define it in a clear and same way for both instruments here. By "frequency vector", we mean a vector containing the frequency value corresponding to each channel of the FFTS.

## Quantities of interest
There are different quantities needed for defining the frequency vector for an instruments and these are:
* The number of channels from the FFTS (N)
* The frequency of the local oscillator(s) used for the downconversion of the signal (fLO1,2,...)
* The sampling rate of the FFTS (FS)

## Frequency Vector
From the above mentionned quantities we should be able to derive the frequency vector corresponding to the FFTS channels. 

### Bandwidth
The measured range for a FFTS is [0 FS/2]. In the case of the Aquiris 240, FS= 2GHz giving a bandwidth of 1GHz.

### Frequency resolution
The frequency resolution (df) is then equal to the measured frequency range divided by the number of bins => df=(FS/2)/N 

In our case we get:
* SOMORA: 61.035 kHz
* GROMOS: 30.518 kHz

### Frequency axis
To compute the frequency axis we then use the fact that the first channel corresponds to 0Hz (DC channel) and the last one to (FS/2)*(1-1/N):

if=(FS/2)*[0:1/N:(1-1/N)]

### "absolute frequency"
To find the absolute frequency axis, we have to make use of the local oscillators frequencies.

#### SOMORA:
It has 3 local oscillators:
1. fLO1 = 149.275GHz (2*74.6375) 
1. fLO2 = 5.6GHz
1. fLO3 = 2GHz

It is LSB:
DC channel corresponds to the frequency: fLO1 - fLO2 - fLO3 = 141.675GHz

#### GROMOS:
It has 2 local oscillators:
1. fLO1 = 145.875GHz 
1. fLO2 = 3.6GHz

TO CHECK !!!!


It is LSB:
The second local oscillator is downconverting to [-500MHz;500MHz].
DC channel corresponds to the frequency: fLO1 - fLO2 - 500MHz = 141.775GHz
