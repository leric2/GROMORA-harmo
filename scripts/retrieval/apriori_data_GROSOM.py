#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 

@author: eric

Collection of functions for dealing with GROSOM a-priori data

Using the Atmosphere Class to define a 

Including : 
    * a-priori data
"""

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

from retrievals import arts
from retrievals.data.ecmwf import ECMWFLocationFileStore
from retrievals.data.ecmwf import levels

class APrioriDataGROSOM(arts.Atmosphere):
    '''
    TO DO ?
    '''
    def __init__(self, apriori_type):
        self.type = apriori_type

        pass

def extract_ecmwf_ds(ECMWF_store_path, t1, t2):
    '''
    Building the ecmwf store for atmospheric state
    '''
    ecmwf_store = ECMWFLocationFileStore('/home/eric/Documents/PhD/GROSOM/ECMWF', 'ecmwf_oper_v2BERN_%Y%m%d.nc')
    ecmwf_ds = (
        ecmwf_store.select_time(t1, t2 )
        .mean(dim='time')
        .swap_dims({"level":"pressure"})
    )
    return ecmwf_ds

def get_apriori_fascod(retrieval_param):
    fascod_atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])
    return fascod_atm

def get_apriori_atmosphere(retrieval_param):
    '''
    Defining a-priori atmosphere from ECMWF operationnal dataset

    First, we initiate an Atmosphere from Fascod climatology (not needed) before
    updating the required fields from ECMWF operationnal dataset.

    TO DO: add the CIRA86 cimatology, otherwise the ecmwf data do not go high enough.
    '''
    
    ECMWF_store_path = retrieval_param['ecmwf_store_path']

    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, t1 = '2019-04-16 03:30', t2 = '2019-04-16 16:45')

    atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])

    # Temperature
    atm.set_t_field(ds_ecmwf['pressure'].values, ds_ecmwf['temperature'].values)

    # z field
    atm.set_z_field(ds_ecmwf["pressure"].values, ds_ecmwf["pressure"].values)

    # Ozone
    atm.set_vmr_field(
        "O3", ds_ecmwf["pressure"].values, ds_ecmwf['ozone_mass_mixing_ratio'].values
        )
    
    return atm