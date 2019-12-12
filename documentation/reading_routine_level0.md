# Reading routine

## Objective and role of this function
Reading raw data from the instrument for a given period of time and store it in standardized matlab structure.

If one can not store the complete day in 1 structure, we will have to include quality check of the raw data directly in this scripts.

### Calls from
Main routine

### Calling
None

### Previously
Variable definition

### Next
Quality control of the raw data

## Inputs
* Location of the raw data
* Time period (1 day ?)
* Instruments name ?

## Outputs
* 1-N standardized matlab structure containing the **raw** data from 1-N instruments.
* 1-N standardized matlab structure containing the **log** data.
* Error code if the raw data are not found (or only partly) for this day.

## Structure
For each day, 2 files to read for each instrument:

* 1 .txt daily file with the following structure:
    * Comment lines starting with %
    * 1 line with the name of the  M parameters (M is different for each instrument)
    * Arbitrary number of line N (1 per timestamp)
    
* 1 .bin daily file:
    * 32bit floating point data of size (N x #channels)

## How it was done in the past
### GROMOS
#### Log file 
Read and saved directly as .mat file

```Matlab
y = readtext([file '.txt'], '\t');
[M,N] = size(y(2:end,:));    % M = number of spectra, N = number of housekeeping data
```
#### Binary 
Read line by line. Select first the lines depending on their indices (h,c,a) before reading the spectrum to perform directly the calibration.

### SOMORA
#### Log file 
Same as GROMOS

#### Binary 
Also goes line by line and calibrate it direclty

## Questions remaining
Is it possible to keep the full day for both instruments before going to the next step ? 

It is done so by Axel so it might by possible to apply this ? --> to test on both instruments

```Matlab
% read complete binary data in one array (huge!)
if nargout>1
    fid = fopen( [file '.bin'], 'r', 'b');
    data = fread(fid, [channels,M], 'float32');
    fclose(fid);
end
```
