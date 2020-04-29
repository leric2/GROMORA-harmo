#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Module reading level1b data 

Example:
    E...

        $ python example_google.py

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * everything

"""

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

def read_level1b(filenameLevel1b):
    """Example function with types documented in the docstring.
    Description HERE

    Args:
        param1 (int): The first parameter.
        param2 (str): The second parameter.

    Returns:
        bool: The return value. True for success, False otherwise.

    """
    
    DS = xr.open_dataset(
        filenameLevel1b + ".nc",
        group="spectrometer1",
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        #use_cftime=True,
        )

    globalAttributes=xr.open_dataset(filenameLevel1b + ".nc").attrs
    
    METEO=xr.open_dataset(
        filenameLevel1b+".nc",
        group="meteo",
        decode_times=True,
        decode_coords=True,
        )
    
    return DS,METEO,globalAttributes

def main():
    filename="/home/eric/Documents/PhD/GROSOM/Level1/SOMORA_level1b_AC240_2019_04_16"
    level1b,METEO,globalAttributes=read_level1b(filename)
    return level1b

if __name__=="__main__":
    main()