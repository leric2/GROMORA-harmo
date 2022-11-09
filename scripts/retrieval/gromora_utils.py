#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Utilities

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

Todo: all

"""
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

K_B = 1.38065e-23 #[J/K]
MAIR = 28.9644 #[g/mol]
MO3 = 47.9982 #[g/mol]
DU= 2.69e20 #[molec/m2]
R = 287.058 #[J/kgK]

def save_single_pdf(filename, figures):
    """
    Save all `figures` to a single PDF. taken from Jonas
    """
    with PdfPages(filename) as pdf:
        for fig in figures:
            pdf.savefig(fig)

def save_pngs(basename, figures):
    """
    Save all `figures` to a single PDF. taken from Jonas
    """
    for i, fig in enumerate(figures):
        filename = basename+str(i)+'.png'
        fig.savefig(filename, dpi=600, facecolor='w', edgecolor='w')


def var_allan(y, axis=0):
    """
    Compute Allan variance of `y` along `axis`.
    """
    var = np.mean(np.square(np.diff(y, axis=axis)), axis=axis) / 2
    return var


def sideband_response_theory(RF, delta_z, polarisation_change=True):
    """ Simulation of the theoretical sideband response from a MPI

    Args:
        RF (numpy array): frequency vector, in Hz
        delta_z (double): _description_
        polarisation_change (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """    
    if polarisation_change:
        return np.sin(np.pi*delta_z*RF/3e8)**2
    else:
        return np.cos(np.pi*delta_z*RF/3e8)**2

def o3_vmr2mmr(o3_vmr):
    """function to convert volume to mass mixing ratios of ozone

    Args:
        o3_vmr (_type_): ozone volume mixing ratio

    Returns:
        o3_mmr: ozone mass mixing ratio
    """       
    o3_mmr = o3_vmr*MO3/MAIR
    return o3_mmr

def density(p, T):
    return p/(R*T)

def o3_mmr2mass_conc(o3_mmr,p,T):
    return o3_mmr*density(p, T)

def o3_vmr2number_density(o3_vmr,p,T):
    return o3_vmr*p/(K_B*T)