# Calibration routine

## Objective and role of this function

### Calls from


### Calling


### Previously
...

### Next
...

## Inputs


## Outputs
* Calibrated spectra in a structure



## Structure

Reformat the spectrum

Optional

Transform the raw data (vector) into a matrix.

Level 0 checks

Optional

Perform a few checks on the raw data. Saved as an attributes in the level 1a.

\paragraph{Flip the spectrum}

Conditional

For some instrument, flipped spectrum

Inversion

\paragraph{Plot the raw counts}

Optional

Simple function plotting (uglily) the raw counts. The three additional inputs (lowerLim, upperLim, N) are just defining
the lower/upper limit of the FFTS counts and the number (N) of raw spectra to plot distributed regularly on the whole day.