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
        RF (_type_): _description_
        delta_z (_type_): _description_
        polarisation_change (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """    
    if polarisation_change:
        return np.sin(np.pi*delta_z*RF/3e8)**2
    else:
        return np.cos(np.pi*delta_z*RF/3e8)**2

