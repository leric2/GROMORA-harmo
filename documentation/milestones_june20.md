# GROSOM Milestones

## Calibration
* Discuss with Franzisca about MIAWARA implementation and feedbacks
* Update the calibration routine for GROSOM and MOPI
* Run some raw data processing (1 month ?) to spot missing stuff

## Retrieval
* Check CPU usage by ARTS...
* Use pyarts instead of Python in pyretrievals
* Sensors definition and COVMAT error...
* Spectroscopy...! 
* Think about the structure and outputs, 1 file per day ? per hour ?

## Atmospheric data
* Download ERA5 data ? Global files ?
* Pre-processing of ECMWF data with pyretrievals -> to do with netCDF oper V2 !
* Store the apriori datasets separated from the level2 files and read it when needed -> in the form of the ECWF store ?
* Include H20 a-priori data for the retrieval of water vapor --> merge ECMWF with Fascod ?
* Decide on which Ozone a-priori data --> some literature reading here ?

## Presentation ideas
* The calibration routine --> make some nice plots
* Data structure for level1a and level1b
* Retrievals organisation
* pyretrievals tools, typhon, pyarts and the future ?
* ARTS3 ?
* Workflow for the GROSOM projet, how to deal with the matlab/python interfaces
* first retrievals 