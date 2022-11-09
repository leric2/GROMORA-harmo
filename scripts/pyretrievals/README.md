# pyretrievals

Building on the great projects [ARTS](http://radiativetransfer.org), 
[Typhon](https://radiativetransfer.org/misc/typhon/doc/index.html) and now relying on [PyARTS](https://pypi.org/project/pyarts/)
the pyretrieval library provides some tools written by Jonas Hagen to run simulations and retrievals for the WIRA-C instrument.

The initial library is available on [Github](https://github.com/jonas-hagen/pyretrieval) and was written for a development version of ARTS. 
Since ARTS has now been updated and Jonas is now gone, the pyretrievals library has been updated and integrated directly within the 
GROMORA project repository so that it is easier for non-Python specialist to use it.

## Environment setup

Using conda (recommended), you can build your environment using the environment file provided in scripts: [env_file_GROMORA.txt](../env_file_GROMORA.txt)

To use ARTS, you also need to specify some paths in your environment. 
The variables `$ARTS_BUILD_PATH` and `$ARTS_DATA_PATH` `ARTS_INCLUDE_PATH` should be either exported in your shell or loaded at the beginning of your scripts.

One way to do it is with the use of dotenv package in python:
```
load_dotenv('/opt/anaconda/.env.birg-arts24')

ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']
```

Example for the shell export:

```bash
ARTS_DATA_PATH=/opt/arts-dev/arts-xml-data/
ARTS_SRC_PATH=/opt/arts-dev/arts/
```

## Examples and tests

The ``examples`` directory contains a ``test_oem.py`` file which is similar to the ARTS cfile
``controlfiles/artscomponents/test_oem.arts``. It simulates the ozone line at 110 GHz and retrieves ozone VMR, frequency shift and a polynomial baseline.

This should be used as a test to see if the package works fine.
