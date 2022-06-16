# retrieval_param

## Summary
This is a Python dictionary containing all the information to launch a retrievals for
a GROMORA intrument. It contains all parameters required for the retrievals of a specific instrument and day.

---

## Table of Contents
1. [Building calibrationTool](#building-calibrationtool)
2. [Parameters](#parameters)

---

## Building retrieval_param

The *retrieval_param* dictionary is unique for each instrument and each date of processing. 

The building of the *retrieval_param* structure is done at different places
within the GROMORA retrieval routine:

### 1. *gromora_retrievals*

The creation of the *retrieval_param* is mainly defined within the *define_retrieval_param* function located in the **DataRetrieval** class from the *gromora_retrievals*.

### 2. instrument classes

All the instrument specific parameters (spectrometer, sideband frequencies, ...) are defined with the instrument classes.

### 3. Main script

When *retrieval_param* is initialized, some important parameters are defined directly within the main script of the GROMORA retrieval (*retrieve_all.py*). These are generally mostly boolean parameters or key parameters that are defined there so that they can be changed easily.

### 4. During the execution

After the first 3 steps, the building of the *retrieval_param* structure is finished.

Only a few specific parameters will be added during the processing of the routine as they will be read from the level 1 file (e.g. zenith angle).

---

## Parameters

### General parameter for ARTS

|variable | type | Description |
|-------|------|:-----------|
|ARTS_DATA_PATH| str | path for ARTS| 
|ARTS_BUILD_PATH|str | path for ARTS |
|ARTS_INCLUDE_PATH|str | path for ARTS|

---

### Main parameters

|variable | type | Description |
|-------|------|:-----------|
|retrieval_type| int | type of retrieval to do 
|retrieval_quantities| str |o3_h2o_fshift_polyfit_sinefit|
| obs_freq | double | observation frequency in Hz 
|date | Timestamp | the day of the retrieval
| integration_cycle | int | the integration cycle to retrieve 
| time| numpy array of datetime64 | the weighted time of measurement
| time_start | numpy array of datetime64 | start time of measurement 
| time_stop| numpy array of datetime64 | stop time of measurement
| oem_method | str | type of OEM algorithm 
| max_iter| double | max number of iteration
| stop_dx|  double | threshold to reach for dx
| lm_ga_setting| list of double | parameter to use in case we use the LM method for the OEM

---

### Plotting and saving variables
| variable | type  | Description |
|---------|------|:-----------:|
|FM_only| boolean | 1 to only do the forward model
|show_FM| boolean| 1 to show the FM
|verbose | int | verbosity of the printing 
|plot_meteo_ds| boolean |  1 to plot the meteo variables | 
|show_f_grid| boolean | 1 to show the frequency simulation grid
|plot_opacities| boolean | 1 to plot the Ingold opacities
|plot_o3_apriori_covariance| boolean | 1 to plot the ozone a priori covariance matrix

---

### Geolocation data 

|variable | type | unit | description | 
|------|------|------|:-----------|
| lon | array | degree_north | latitude defined according to WGS84
| lat | array | degree_east | longitude defined according to WGS84
| surface_altitude | double | meter | altitude of the surface above see level 
| observation_altitude | double | meter | altitude of the observation above see level 
| station_altitude | double | meter | real altitude of the station above see level 
| zenith_angle | double | degree measured clockwise positive, 0 deg is zenith | zenith observation angle 
| azimuth_angle | double | degree measured clockwise positive, 0 deg is northwise | azimuth observation angle 
| pointing_angle_corr | double | degree  | potential offset of the zenith angle


---

### Atmosphere and a priori

| variable | type  | Description |
|---------|------|:----------|
| atm | str | type of atmosphere to use for the retrieval
| ptz_merge_method | str | method for merging ECMWF and CIRA86  
| ptz_merge_max_Tdiff | double | max temperature difference to accept for merging ECMWF and CIRA86 profiles
| ecmwf_store_location | str |  path to ECMWF files
| ecmwf_prefix | str | prefix of ECMWF files
| cira86_path | str | path to CIRA montly climatology folder
| extra_time_ecmwf | double | extra time to consider for ECMWF profile (in hour)
| selected_species | list of str | list of all atmospheric species to consider in the forward model  
| o3_apriori | str | type of a priori to take for ozone
| o3_apriori_covariance | str | type of ozone covariance to consider 
| water_vapor_model | str | absorption model for water vapor  | 
| h2o_apriori | str | type of a priori to take for water vapor
| apriori_O3_cov| double | variance for ozone
| apriori_H2O_stdDev| double  |relative std deviation for water vapor (in %) | 
| o2_model | str | absorption model for oxygen | 
| n2_model | str | absorption model for nitrogen | 
| spectroscopy_type | str | type of spectoscopy | 
|line_file | str | the spectroscopy file to use for the retrievals  
| waccm_file | str | filename for WACCM climatology
|p_surface | array | surface pressure (hPa)
|T_surface | array | surface temperature (K)

---

### Simulation and retrieval grids

| variable | type  | Description |
|---------|------|:----------|
| retrieval_grid_type | str | type of retrieval grid
| z_top_sim_grid | double | altitude max for FM (in meter)
| z_bottom_sim_grid| double | altitude min for FM (in meter)
| z_resolution_sim_grid| double | altitude resolution for FM (in meter)
| z_top_ret_grid| double| altitude max for retrieval (in meter)
| z_bottom_ret_grid| double| altitude min for retrieval (in meter)
| z_resolution_ret_grid| double| altitude max resolution retrieval (in meter)
| retrieval_h2o_grid_type | str | type of water vapor retrieval grid
| h2o_pressure| numpy array | retrieval grid for water vapor continuum

---

### Sensors parameters

| variable | type  | Description |
|---------|------|:----------|
| sensor | str | type of sensor to simulate
| SB_bias| double | bias on the sideband
| sideband_response | str | way to simulate the sideband response
| use_all_channels| boolean | 1 to use all FFT channels in the measurement vector
| window_corrected_spectrum|  boolean  | 1 to use the window corrected spectra
| f_max |  double | maximum frequency of the measurement vector (Hz)
| f_min |  double | minimum frequency of the measurement vector (Hz)
| f_shift| double | offset on the frequency vector (in Hz) 
| number_of_freq_points| double | number of simulation frequency points
| irregularity_f_grid| double | irregularity factor of the frequency grid (refined at line center)
| bandwidth| double | bandwidth of the spectrometer in Hz 
| noise_covariance | str | method to use for defining the noise covariance matrix
| increased_var_factor| double | factor by which to multiply the noise variance 
| binned_ch| boolean | 1 to bin the measurement spectrum (not tested)

---


### Instrumental baseline retrievals

| variable | type  | Description |
|---------|------|:----------|
| poly_order| double | degree of polynomial to retrieve 
| covmat_polyfit_0 | double | variance for polynomial retrieval degree 0
| covmat_polyfit_1 | double | variance for polynomial retrieval degree 1
| covmat_polyfit_2 | double | variance for polynomial retrieval degree 2
| sinefit_periods| numpy array | periods for sine baselines retrievals 
| sinefit_covmat| numpy array | variance to use for sine baselines retrievals 





