# Structure level 1a data
For detailed information, see https://www.unidata.ucar.edu/software/netcdf/docs/user_guide.html
## Format
It was agreed to use the netCDF format for storing all levels of our data. There are 4 different types of base format for netCDF which can be divided between the ***netCDF classic base format*** and the enhanced ***netCDF-4/HDF5***.

The ***netCDF classic base format*** include CDF-1, CDF-2 and CDF-5. From version 4.0, the ***netCDF-4/HDF5*** is using HDF5 format as a storage layer (which make it readable by HDF5) and offers different advantages compared to the older classical format in terms of dimensions, group definitions, etc... This is the base format that we will use for storing the level1a and the rest of the relevant data. 

## Data Model
netCDF classic dataset is stored as a single file containing 2 parts:
* header: information about dimensions, attributes and variables which all have both a name and an ID by which they are identified.
* data: fixe- or variable-size data for variables that have (un)limited dimension

### Dimensions
It represents either a real physical quantity or index other quantities. It has both a name and a length (positive int., can be unlimited)

### Variables
Array of values of the same type. It has a name, a data type and a shape (described by its list of dimensions).

It can have attributes.

#### Coordinate variables
These are one-dimensional variable with the same name as its dimension. 

### Attributes
These are the metadata and store information about the data. They can be either ***global*** and store information about the dataset or ***local*** for a specific variable.  

"Attributes are more dynamic than variables or dimensions; they can be deleted and have their type, length, and values changed after they are created, whereas the netCDF interface provides no way to delete a variable or to change its type or shape."

* Global attributes are identified by its name and a special "global variable"
* Local attributes are identified by its name and the name (or ID) of the specific variable

## Conventions
There are multiple conventions and good practice for writing netCDF files. For the sake of compatibility, we will start from existing conventions to write our netCDF file (see http://cfconventions.org/) and we will adapt it if needed.