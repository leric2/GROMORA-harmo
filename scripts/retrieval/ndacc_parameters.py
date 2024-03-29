#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 

@author: eric

Just a file storing the NDACC parameters needed to write the HDF4 files for GROMORA.
"""
import sys, os, time
from os.path import dirname, abspath, join
import numpy as np

SDS_Name_List = ['LATITUDE.INSTRUMENT',
                 'LONGITUDE.INSTRUMENT',
                 'ALTITUDE.INSTRUMENT',
                 'DATETIME',
                 'ANGLE.VIEW_AZIMUTH',
                 'ANGLE.VIEW_ZENITH_MEAN',
                 'ANGLE.SOLAR_ZENITH_MEAN',
                 'OPACITY.ATMOSPHERIC_EMISSION',
                 'DATETIME.START',
                 'DATETIME.STOP',
                 'INTEGRATION.TIME',
                 'ALTITUDE',
                 'PRESSURE_INDEPENDENT',
                 'TEMPERATURE_INDEPENDENT',
                 'O3.MIXING.RATIO.VOLUME_EMISSION',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.RANDOM.STANDARD',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.SYSTEMATIC.STANDARD',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.COMBINED.STANDARD',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_RESOLUTION.ALTITUDE',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_APRIORI',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_APRIORI.CONTRIBUTION',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_AVK',
                 'O3.COLUMN.PARTIAL_EMISSION',
                 'O3.NUMBER.DENSITY_EMISSION']

SDS_DataType_List = ['REAL',
                     'REAL',
                     'INTEGER',
                     'DOUBLE',
                     'REAL',
                     'REAL',
                     'REAL',
                     'REAL',
                     'DOUBLE',
                     'DOUBLE',
                     'DOUBLE',
                     'INTEGER',
                     'REAL',
                     'REAL',
                     'REAL',
                     'REAL',
                     'REAL',
                     'REAL',
                     'INTEGER',
                     'REAL',
                     'REAL',
                     'REAL',
                     'REAL',
                     'REAL']

VAR_DESCRIPTION_List = ['Latitude of the observation site (deg)',
                        'Longitude of the observation site (deg)',
                        'Altitude of the observation site (m)',
                        'Representative Date and Time of the Measurement (MJD2K)',
                        'Viewing Azimuth angle (deg)',
                        'Mean Viewing Zenith angle (deg)',
                        'Mean Solar Zenith angle of measurement (deg)',
                        'Tropospheric Optical Depth at 142 GHz (Np)',
                        'Start Date and Time of the Measurement (MJD2K)',
                        'Stop Date and Time of the Measurement (MJD2K)',
                        'Effective Measurement Time (h)',
                        'Altitude as geometric height (m)',
                        'Derived Atmospheric Pressure (hPa)',
                        'Derived Temperature Profile (K)',
                        'O3 Volume Mixing Ratio Profile (ppmv)',
                        'O3 Profile Observation Error (ppmv)',
                        'O3 Profile Smoothing Error (ppmv)',
                        'O3 Profile Total Error (ppmv)',
                        'O3 Mixing Ratio Altitude Resolution (m)',
                        'A priori O3 Profile (ppmv)',
                        'Contribution from the A priori to the retrieved O3 Profile (%)',
                        'O3 Retrieval Averaging Kernels',
                        'O3 partial column density (DU)',
                        'O3 number density (molecules m-3)'];
                        
VAR_NOTES_List = [' ',
                  ' ',
                  ' ',
                  'Averaged time of sky measurements',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  'Merge of ECMWF operational analysis and CIRA86',
                  'Merge of ECMWF and CIRA86',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  'Resolution defined as FWHM of the averaging kernel peak for that altitude level.',
                  'Monthly day and night profile extracted from WACCM free running model',
                  ' ',
                  ' ',
                  ' ',
                  ' ']
                    
VAR_DEPEND_List = ['CONSTANT',             
                   'CONSTANT',
                   'CONSTANT',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE;ALTITUDE',
                   'DATETIME',
                   'DATETIME;ALTITUDE']
               
VAR_UNITS_List = ['deg',
                  'deg',
                  'm',
                  'MJD2K',
                  'deg',
                  'deg',
                  'deg',
                  'Np',
                  'MJD2K',
                  'MJD2K',
                  'h',
                  'm',
                  'hPa',
                  'K',
                  'ppmv',
                  'ppmv',
                  'ppmv',
                  'ppmv',
                  'm',
                  'ppmv',
                  '%',
                  '1',
                  'DU',
                  'molec m-3']
    
VAR_SI_CONVERSION_List = ['0.0;1.74533E-2;rad',
                          '0.0;1.74533E-2;rad',
                          '0;1;m',
                          '0.0;86400.0;s',
                          '0.0;1.74533E-2;rad',
                          '0.0;1.74533E-2;rad',
                          '0.0;1.74533E-2;rad',
                          '0.0;1.0;1',
                          '0.0;86400.0;s',
                          '0.0;86400.0;s',
                          '0.0;3600.0;s',
                          '0;1;m',
                          '0.0;1.0E2;kg m-1 s-2',
                          '0.0;1.0;K',
                          '0.0;1.0E-6;1',
                          '0.0;1.0E-6;1',
                          '0.0;1.0E-6;1',
                          '0.0;1.0E-6;1',
                          '0;1;m',
                          '0.0;1.0E-6;1',
                          '0.0;0.01;1',
                          '0.0;1.0;1',
                          '0.0;4.4614E-4;mol m-2',
                          '0.0;1.66054E-24;mol m-3']

VAR_VALID_MIN_List = [-90.0,
                      -180.,
                      0.,
                      np.nan,
                      0.,
                      0.,
                      0.,
                      0.,
                      np.nan,
                      np.nan,
                      0.,
                      0,
                      0.,
                      100.,
                      -0.1,
                      -0.1,
                      -0.1,
                      -0.1,
                      1000,
                      -0.1,
                      -70.,
                      -1.,
                      0.,
                      0.]

VAR_VALID_MAX_List = [90.,
                      180.,
                      5000.,
                      np.nan,
                      360.,
                      180.,
                      180.,
                      5.,
                      np.nan,
                      np.nan,
                      1000.,
                      120000,
                      1100.,
                      500.,
                      20.,
                      20.,
                      20.,
                      20.,
                      35000,
                      20.,
                      200.,
                      1.,
                      500.,
                      1.e19]

list_indices = {
    'DATETIME':3, 
    'ALTITUDE':11, 
    'LATITUDE':0,
    'LONGITUDE':1, 
    'ALTITUDE_INSTRUMENT':2,
    'ANGLE_VIEW_AZIMUTH':4,
    'ANGLE_VIEW_ZENITH':5,
    'ANGLE_SOLAR_ZENITH_MEAN':6,
    'OPACITY_ATMOSPHERIC_EMISSION':7,
    'DATETIME_START':8,
    'DATETIME_STOP':9,
    'INTEGRATION_TIME':10,
    'PRESSURE_INDEPENDENT':12,
    'TEMPERATURE_INDEPENDENT':13,
    'O3_MIXING_RATIO_VOLUME_EMISSION':14,
    'O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_RANDOM_STANDARD':15,
    'O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_SYSTEMATIC_STANDARD':16,
    'O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_COMBINED_STANDARD':17,
    'O3_MIXING_RATIO_VOLUME_EMISSION_RESOLUTION_ALTITUDE':18,
    'O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI':19,
    'O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI_CONTRIBUTION':20,
    'O3_MIXING_RATIO_VOLUME_EMISSION_AVK':21,
    'O3_COLUMN_PARTIAL_EMISSION':22,
    'O3_NUMBER_DENSITY_EMISSION':23
    }
