# Level 1

## Table of Contents
1. [Format](#netcdf-format)
1. [Tools and Conventions](#tools-and-conventions)
3. [File structure](#grosom-file-structure)

    1. [Global attributes](#global-attributes)
    2. [Group 1: spectrometer](#spectrometer1-group)
    2. [Group 2: flags](#flags-group)
    2. [Group 3: meteo](#meteo-group)

## NetCDF Format

The Network Common Data Form (netCDF) is a data format that allows creation,
access and sharing of scientific datasets. It is a standard within the
scientific communit for the storage and exchange of array-oriented data which
means that there are multiple tools to deal with this format (Panoply, ncdump,
...) and conventions on how to write netCDF files.

There are multiple conventions and good practice for writing netCDF files. For
the sake of compatibility, we will start from existing conventions to write our
netCDF file (see \href{http://cfconventions.org/}{CF conventions}) and we will
adapt it if needed.

The advantages of the netCDF format is that it enables to store and access data
in an efficient way and which is machine-independant. Moreover, many API exist
for writing and reading netCDF files and it is especially 


It was agreed to use the netCDF format for storing all levels of our data. There are 4 different types of base format for netCDF which can be divided between the ***netCDF classic base format*** and the enhanced ***netCDF-4/HDF5***.

The ***netCDF classic base format*** include CDF-1, CDF-2 and CDF-5. From version 4.0, the ***netCDF-4/HDF5*** is using HDF5 format as a storage layer (which make it readable by HDF5) and offers different advantages compared to the older classical format in terms of dimensions, group definitions, etc... This is the base format that we will use for storing the level1a and the rest of the relevant data. 

GEOMS standards used by the NDACC : https://gs614-avdc1-pz.gsfc.nasa.gov/index.php?site=2024100160. I think that those are not so nice and that we shoud better stick with the cf conventions (see later).

### Data Model
netCDF classic dataset is stored as a single file containing 2 parts:
* header: information about dimensions, attributes and variables which all have both a name and an ID by which they are identified.
* data: fixe- or variable-size data for variables that have (un)limited dimension

### Dimensions
It represents either a real physical quantity or index other quantities. It has both a name and a length (positive int., can be unlimited)

The level 1 files for the GROSOM project have the time as the main dimension. On 

### Variables
Array of values of the same type. It has a name, a data type and a shape (described by its list of dimensions).

#### Coordinate variables
These are one-dimensional variable with the same name as their dimension. 

### Attributes
These are the metadata and store information about the data. They can be either ***global*** and store information about the dataset or ***local*** for a specific variable.  

"Attributes are more dynamic than variables or dimensions; they can be deleted and have their type, length, and values changed after they are created, whereas the netCDF interface provides no way to delete a variable or to change its type or shape."

* Global attributes are identified by its name and a special "global variable"
* Local attributes are identified by its name and the name (or ID) of the specific variable

For detailed information, see https://www.unidata.ucar.edu/software/netcdf/docs/user_guide.html

## Tools and Conventions
There are multiple conventions and good practice for writing netCDF files. For the sake of compatibility, we will start from existing conventions to write our netCDF file (see http://cfconventions.org/) and we will adapt it if needed.

In terms of tools, they are plenty of possibilities to deal with netCDF data files. Some worth mentionning are the excellent [Panoply](https://www.giss.nasa.gov/tools/panoply/) and the netCDF utilities from Unidata. The latter enable to have a very quick look at the data or to manipulate (copy, append, extract subset) netCDF files very easily.

## GROSOM file structure

GROSOM level 1 files are stored in netCDF-4 format.

Data are stored in daily files with time interval corresponding to the calibration (level 1a) or intergation (level 1b) time. Both levels have the same structure and more less the same variables. All files contains the following elements:
* Global attributes: some meta information about the file content
* Three 3 groups of variables:
    * spectrometer1: contains all outputs from the calibration/integration from the main spectrometer in the instrument.
    * flags: contains a set of flags to assess the quality of the data stored in spectrometer1
    * meteo: contains all meteorological variables at the location of the instrument.

Note that each groups and variables within have in addition its own set of attributes. 

The details of the level 1 for the GROSOM project is presented below:

### Global attributes

#### Level 1a
| Attributes | type  | Description |
|------|------|:-----------|
| title | str  | Content of the file |
| location | str  | location of the instrument |
| source | str  | type of data (see [NDACC](https://gs614-avdc1-pz.gsfc.nasa.gov/index.php?site=1925698559)) |
| name | str  | name of the PI |
| institution | str  | instituion responsible for this instrument |
| contact | str  | address of the institution |
| mail | str  | email address of the PI |
| instrument | str  | name of the instrument |
| number_of_spectrometer | double  | number of spectrometer saved in this file |
| history | str  |  |
| references | str  |  |
| comment | str  | miscellaneous comments |
| raw_filename | str  | file name of the raw data |
| raw_data_software_version | str  |  |
| calibration_version | str  |  |
| raw_file_comment | double  |  |
| raw_file_warning | str  | warning on the raw files |
| labview_logfile_warning | str  | check if a log entry was present for this day in the labview log |
| data_start_date | double  | first datetime in MJD2K of this file |
| data_stop_date | double  | last datetime in MJD2K of this file |
| filename | str  | complete filename of this file |
| creation_date | str  | creation date of this file |
| featureType | str | [CF conventions](http://cfconventions.org/) |

#### Level 1b

For level 1b, the global attributes are moreless the same only with the following additions:

| Attributes | type  | Description |
|------|------|:-----------|
| filename_level1a | str  | name of the level 1a file |
| creation_date_level1a | str  | creation date of the level 1a file |

### spectrometer1 group

The main group for both level 1 files is named spectrometer1. It is where we store all variables extracted during the calibration process. 

As attribute for this group : spectrometer_type

### Dimensions:
1. time (unlimited)
2. channel_idx

### Coordinate variables:
| Coordinates | type | units | other attributes | Description |
|------|------|------|------|:-----------|
| time | double | days since 2000-01-01 00:00:00 | calendar | mean time recorded at the beginning of all sky measurements during this calibration cycle |
| channel_idx | long |- | - | indices of the channels for this spectrometer|

### Variables:

All variables should contains the following attributes ([CF conventions](http://cfconventions.org/)): units, standard_name, long_name, _FillValue and description. 

| variables | type | dimension | units | standard_name | long_name | description |
|------|------|------|------|------|------|:-----------|
| lat | float | time | degree_north | latitude | station latitude | latitude defined according to WGS84 |
| lon | float | time| degree_east | longitude | station longitude | longitude defined according to WGS84 |
| alt | float | time | meter | atitude | station atitude | above see level |
| azimuth_angle | float | time | degree | sensor_azimuth_angle | azimuth angle | angle measured clockwise positive, 0 deg is northwise |
| MJD2K | double | time | MJD2K | - | -| mean time recorded at the beginning of all sky measurements during this calibration cycle |
| year | long | time | - | - | - | - |
| month | long | time | - | - | - | - |
| day | long | time | - | - | - | - |
| time_of_day  | double | time | hour | TOD | time of day | Time of the day |
| first_sky_time  | double | time | days since 2000-01-01 00:00:00 | - | - | time of the first sky measurements in this calibration cycle |
| last_sky_time  | double | time | days since 2000-01-01 00:00:00 | - | - | time of the last sky measurements in this calibration cycle|
| time_min  | double | time | days since 2000-01-01 00:00:00 | - | - | minimum theoretical start time for this calibration cycle" |
| Tb  | double | time, channel_idx | Kelvin | brightness_temperature | Tb | calibrated brightness temperature for this cycle |
| stdTb  | double | time, channel_idx | Kelvin | std_Tb | standard variation of Tb | standard deviation of brightness temperature for this cycle |
| frequencies | double | channel_idx | Hz | frequency | frequency vector | frequency vector for the spectrometer |
| intermediate_freq | double | channel_idx | Hz | intermediate_frequency| intermediate frequency vector | intermediate frequency vector for the spectrometer |
| mean_std_Tb  | double | time | Kelvin | mean_stdTb | mean standard variation of Tb | mean standard deviation of brightness temperature for this cycle (without bad channel) |
| THot  | double | time | Kelvin | hot_load_temperature | THot | Mean temperature of the hot load |
| stdTHot  | double | time | Kelvin | std_hot_load_temperature | stdTHot | standard deviation of the hot load temperature |
| TNoise  | double | time | Kelvin | noise_temperature | mean noise temperature | mean noise receiver temperature |
| stdTNoise  | double | time | Kelvin | std_noise_temperature | stdTNoise | standard deviation of the noise receiver temperature |
| calibration_time  | double | time | second | calibration_time | calibrationTime | Time interval used for calibrating the spectra |
| mean_sky_elevation_angle  | double | time | degree | elevation_angle | mean sky angle |mean elevation angle of the sky observation during this cycle |
| TRoom  | double | time | Kelvin | room_temperature | TRoom | mean room temperature |
| stdTRoom  | double | time | Kelvin | standard_room_temperature | stdTRoom | standard deviation of room temperature |
| TWindow  | double | time | Kelvin | window_temperature | TWindow | mean window temperature |
| TOut  | double | time | Kelvin | outside_temperature | TOut | mean outside temperature |
| noise_level  | double | time | Kelvin | noise_level | std(diff(Tb)) | describes how noisy is the spectra |
| number_of_hot_spectra  | long | time | - | - | - | ... |
| number_of_cold_spectra  | long | time | - | - | - | ... |
| number_of_sky_spectra  | long | time | - | - | - | ... |
| mean_hot_counts  | double | time | - | - | - | ... |

The same variables are used in the level 1b file with some additions:

| variables | type | dimension | units | standard_name | long_name | description |
|------|------|------|------|------|------|:-----------|
| integration_time  | double | time | second | integration_time | integrationTime | Time interval used for integrating the spectra |
| number_calibrated_spectra  | double | time | - | - | - | ... |
| tropospheric_transmittance  | double | time | - | tropospheric_transmittance | tropospheric transmittance | method |
| tropospheric_opacity  | double | time | - | tropospheric_opacity | tropospheric opacity | method |

## flags group

The flags groups contains the set of flags for the calibrated and integrated spectra contained in the spectrometer1 group. It has the same main dimension and coordinate as the spectrometer1 group and has a secondary dimension corresponding to the number of flags defined.

### Dimensions:
1. time (unlimited)
2. flags

### Coordinate variables:
| Coordinates | type | units | other attributes | Description |
|------|------|------|------|:-----------|
| time | double | days since 2000-01-01 00:00:00 | calendar |mean time of the measurements for this cycle |
| flags | long |- | - | ... |

### Attributes:
It contains 2 attributes, one being the description of the flag group organisation, the other being the number of flags. 

### Variables:
It contains a single variable which is a vector of the defined flags. As attributes, it contains the description of the flag vector, listing which element correspond to which flag.

For the level 1a, we have 6 flags (6 element vector per time stamp):

    :errorCode_1 = "sufficientNumberOfIndices";
    :errorCode_2 = "noiseTemperatureOK";
    :errorCode_3 = "LN2SensorsOK";
    :errorCode_4 = "LN2LevelOK";
    :errorCode_5 = "hotLoadOK";
    :errorCode_6 = "pointingAngleOK";

For the level 1b, we have 2 flags(6 element vector per time stamp):

    :errorCode_1 = "sufficientNumberOfAvgSpectra";
    :errorCode_2 = "tropospheric_transmittance_OK";

For more details on the meaning of these flags, see the [quality control](quality_control_calibration.md) specific documentation. 

### Meteo group

It contains all the meteorological data read during the calibration. In the case of level 1a, we keep the original time stamps from the meteo files and therefore, we get a different time coordinates compared to the spectrometer1 group. For the level 1b, this is not the case as we only consider averaged meteo data.

### Dimensions:
1. time (unlimited)

### Coordinate variables:
| Coordinates | type | units | other attributes | Description |
|------|------|------|------|:-----------|
| time | double | days since 2000-01-01 00:00:00 | calendar |mean time of the measurements for this cycle |

### Variables:

It has 4 main variables with the same attributes as for the spectrometer1 group.

| variables | type | dimension | units | standard_name | long_name | description |
|------|------|------|------|------|------|:-----------|
| air_pressure | double | time | hPa | air_pressure | air pressure | air pressure at the station |
| air_temperature | double | time | Kelvin | air_temperature | air temperature | air temperature at the station |
| relative_humidity | double | time | - | relative_humidity | relative humidity | relative humidity of the air at the station |
| precipitation | double | time | mm | precipitation | precipitation | Accumulation of precipitation during the cycle (from gauge ?) |
