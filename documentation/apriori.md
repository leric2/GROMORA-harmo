# Apriori data for GROSOM (and other) retrievals

## Summary
Description of the workflow for collecting, processing and using apriori data in the GROSOM retrievals. 

There has been an effort from Jonas to harmonize the processing of atmospheric apriori data to work in the retrievals from different radiometers (gromos-c, wira, wira-c, ?). This quick note tries to understand and summarize the methods implemented by Jonas in the pyretrievals package ([pyretrievals documentation](http://www.iapmw.unibe.ch/research/projects/pyretrievals/index.html), section: Atmospheric Data) and apply these to the particular case of GROMOS and SOMORA retrievals. 

## From global dataset to local use for retrieving radiometric data
General idea to go from global dataset download to the apriori profile used in a specific retrievals.

1. Download large global datasets to use in the retrievals (ERA, WACCM) and store the ones that are used on a regular basis on the new atmospheric TUB ? Are the concerned datasets already decided ?

2. Extract and convert the desired dataset into smaller and more practical netcdf files for each location to be stored along the global files for each campaign location --> using the "extract_locations_..." type script. 

3. Use the functions from [ECMWFLocationFileStore](http://www.iapmw.unibe.ch/research/projects/pyretrievals/data.html#retrievals.data.ecmwf.ECMWFLocationFileStore) to access the data at the right time period for performing a given retrieval.

From email with Jonas (11.06.20):
*LocationFileStore is one idea, but better would be a git repository with only such functions, well documented for Matlab and Python (same function names, hence functions).*

This would mean that 1 person should dedicate some time to manage the repos to have something clean and functionnal...

## Using apriori data in the GROSOM retrievals with pyretrievals

Ideas for setting up the apriori atmospheric state for GROSOM retrievals. This is taken from my understanding of what have been done in the WIRA/-C case. 

The definition of the apriori state can be defined into 2 main steps:

- Building an "Atmosphere" that contains all the apriori **raw** fields for temperature, altitude and vmr species. This include the followig steps:

	1. Setup an ArtsController object (with basic required parameters)
	
	2. Setup the simulation grids, *f_grid* and *p_grid*
	
	3. Define a basic atmospheric state based on a climatology (for instance the Fascod using [*Atmosphere.from_arts_xml*](http://www.iapmw.unibe.ch/research/projects/pyretrievals/arts.html?highlight=from_arts_xml#retrievals.arts.Atmosphere.from_arts_xml) method). This creates an Atmosphere object containing already temperature fields, corresponding altitute fields and vmr fields for all species present in the scenario but this is only when we set_atmosphere that the vmr corresponding to the selected absorption species are selected and interpolated with the simulation grids.
	
	4. Extract (for instance ecmwf) apriori profile by using the previously created ECMWFLocationFileStore. 
	
		* Overwrite the temperature and altitude fields in the "basic atmosphere" using the apriori data read from ECMWF or any other apriori.
	
		* Overwrite some selected vmr_fields in the basic atmosphere, typically for the defined absorption species (O3 and H2O in my case). Warning, this steps also applies to the species defined for continuum absorption, meaning N2, O2, ... In the old GROMOS retrieval, the standards Fascod was used for these. 
	
		* Add apriori winds at this point if needed. 

	4. Alternatively, we can set all the required fields in an external apriori dataset and use [*Atmosphere.from_dataset*](http://www.iapmw.unibe.ch/research/projects/pyretrievals/arts.html?highlight=from_dataset#retrievals.arts.Atmosphere.from_dataset) function to convert it into an Atmosphere object. --> We could think about creating a dataset of apriori data for GROSOM and load it at this stage with some thing like GROSOMAprioriDataStore ? In that case, we have to make sure to include all abs_species in our atmosphere vmr_fields raw fields otherwise this will fail when the set_atmosphere is used.
	
- Interpolate the apriori **raw fields** to the simulation grids defined by *p_grid*. This is done using the [*ArtsController.set_atmosphere*](http://www.iapmw.unibe.ch/research/projects/pyretrievals/_modules/retrievals/arts/interface.html#ArtsController.set_atmosphere) function from the ArtsController class.

## Apriori and spectroscopy

For each of the abs_species defined, respectively for each raw_fields included in the apriori data (let's say we want to include also CO2 with O3), we need to set/read some spectroscopic parameters from a file specific to each abs_species ? For WIRA, only O3 spectroscopic parameters were taken into account right ?

*Yes, and only one line! There are thousends of lines, but
selecting the relevant ones and creating a reduced set of lines
enhances performance drastically. This is the approach for all MW
radiometers I guess. No one would ever use the full catalogue (I
guess.)*

To be discussed for GROSOM

## Some questions remaining at that stage
...

