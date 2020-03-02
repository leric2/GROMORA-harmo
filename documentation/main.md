# Main script
## Level0 - level1a (Calibration)
### Summary
Here are some basic information about the routine I have written for GROMOS and SOMORA calibration.

In fact, there are 2 main objects to understand the way I designed this routine:

#### retrievalTool: 
Matlab structure containing every information about the calibration (and then retrieval ??) that we want to run. In the main script (retrieval.m), we are building this structure with 2 types of informations: 

Default or basic parameters: 

The parameters of functions that are inherent to each instrument, for instance the indice that indicates a hot spectra in the log file (ind=2), the number of channel, the localisation of the instruments, ... All those are defined into a function (import_default_retrievalTool.m) that is called in the main script and is building the retrievalTool structure taking as input the instrument name (GROMOS or SOMORA)

"Specific":

Parameters dependant from the type of operation we want to perform, for instance if we want to do the calibration, only the retrieval or if we want to plot some calibrated spectra, ... These have to be defined manually in the main script. 

The idea behind this retrievalTool is that it contains every piece of information for running a successfull retrieval, whatever the instrument. This is idealistic but it enables then to have a standard script following for doing effectively the calibration and the retrieval. Indeed, the process for doing the calibration and the retrieval are then "sequential" and always the same steps if all needed parameters are given. 

This also reduces the number of place where there are instuments dependant properties to only one: the retrievalTool structure. 

Matlab also enables to store functions inside a structure and this enables for instance to use different functions for the same step done on different instruments when it is needed (typically the harmonization of the log information). It enables to avoid the use of ***switch case*** inside a specific functions which would calls for a modification of every ***switch case*** the day we want to integrate new instruments to this routine.

#### run_retrieval.m:
The second main object of this routine is the function doing effectively the calibration. It has only 1 input; the retrivalTool structure and is called from the main script for a given day.

In this function, all steps necessary for a calibration (and then a retrieval) are defined and done sequentially depending on the parameters or functions defined in retrievalTool. 

For the calibration (level0 - level1a), this is executing the following steps:
1. checking the retrievalTool structure
2. reading the raw data
3. harmonizing the log file
4. checking the raw data
5. reformating the spectra
6. flipping the spectra
7. plotting the raw spectra
8. calibrating
9. check and add meta-information to the calibrated spectra
10. plot the calibrated spectra
11. save the calibrated spectra in netCDF daily file

### Error and warning management
Standard way of dealing with error with is to use warnings and/or flags for "unexpected conditions" and errors for "fatal problems" which will stop the scripts.

This part is opened to discussion but what I did until now is:
1. Get some warnings about the raw data, for instance if the size of the binary does not correspond to its log file. --> Those are stored in the netCDF as attributes for each daily files.
2. Some errors, for instance if there is problem when reading the raw data, that completely stops the retrieval for the day.
3. Some flags, linked with every calibration cycle, that are stored along the level1a file in a "vector of flags" with each element corresponding to a specific type of flag (ex: large TSys,...).

### Questions about calibration (to discuss tomorrow)
1. Is calibration routine OK (method) ? Is there a way to check it ?

2. Is the outlier detection OK ? For now, we remove all kind of spectra that have a wrong pointing angle (thresholds ?). For hot and cold spectra, we also check that these are not containing too many channels beyond the median +- n*sigma for each calibration cycle. 

2. Possible improvements: Include tipping curve calibration (approx. every 30 minutes) as additionnal information ? Make some proper error propagation in the calibration formula ? Others ?

3. Thot and Tcold to decide ? (see level1a.md)

4. What other quality flags do we want to add for the calibration cycle ?

5. Output parameters to save as netCDF for the level1a ? See an example of netCDF structure [here](/documentation/example_netCDF/GROSOM_level1a_daily_example.txt).

### Future perspectives

In my opinion, there are different ways we can take for the following steps of the projects:
1. Improve the calibration routine, for instance by including the tipping curves calibration.
2. Prepare for the calibration of the whole time series for both instrument
3. Go to level 1b data
    * Select good calibrated spectra ? Based on what ? flags ? Tb values ?
    * Channel outliers: boxcar filter, absolute value of Tb, value of transmittance ? --> replaced with interpolation ?
    * Windows corrections, others ?
    * Integration of the calibrated spectra on 1h (mean ?)
    * save it in netCDF (at least for debug purposes)
    * continue with Matlab until level 1b (see ARTS mailing list)

<!---
### Suggestion (old)
Take a sort of "object oriented" approach for launching the retrievals. 

We have an object (let's call it retrievalTool) that will contain all the relevant information for the retrieval we want to perform (date, name of the instruments, etc...). We manually define this object for the operation we want to perform (for instance, for running a retrieval between 2 given dates in Payerne, without quality checking anything) and then launch the retrieval through a "run" function.

As launching a retrieval is a step-by-step "linear" operation, we do not need to edit the script effectively doing the step-by-step operation but we only need to edit the retrievalTool object. 

On of the nice thing is that we can store function in this retrievalTool object, so that the retrievals can use different functions for different instruments in the case where it would be needed. This would be defined and documented in the retrievalTool structure and implemented in the main "run" function.

CAUTION: switch cases inside functions should be avoided because it would be painfull to edit them all when a new instruments is added (or a new spectrometers for instance).-->