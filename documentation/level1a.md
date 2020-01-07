# Structure level 1a data
For detailed information, see https://www.unidata.ucar.edu/software/netcdf/docs/user_guide.html
## Format
It was agreed to use the netCDF format. There are 4 different types of netCDF files which can be divided between the ***netCDF classic data format*** and the ***netCDF-4/HDF5***.

The ***netCDF classic data format*** include:
* CDF-1
* CDF-2
* CDF-5

netCDF-4/HDF5: using HDF5 as its base format and compared to the classic netCDF can:
* Enable to use groups and user-defined types
* Can have multiple unlimited dimensions

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
...

### Attributes
These are the metadata and store information about the data. They can be either ***global*** and store information about the dataset or ***local*** for a specific variable.  

"Attributes are more dynamic than variables or dimensions; they can be deleted and have their type, length, and values changed after they are created, whereas the netCDF interface provides no way to delete a variable or to change its type or shape."

* Global attributes are identified by its name and a special "global variable"
* Local attributes are identtified by its name and the name (or ID) of the specific variable



