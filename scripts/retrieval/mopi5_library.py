#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:58:02 2020

@author: eric

Module with function dealing with mopi5 data at every level

"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages


def return_bad_channels_mopi5(number_of_channel, date, spectro):
    '''
    common function for the 2 classes... maybe a better way to do it.
    '''
    if spectro == 'USRP-A':
        bad_channels = np.hstack((np.arange(0,1023), 8092, np.arange(number_of_channel-1023,number_of_channel)))
    elif spectro == 'U5303':
        #bad_channels = np.where(intermediate_frequency > 1000)
        bad_channels = np.hstack((np.arange(0,63), np.arange(11000,number_of_channel)))
    elif spectro == 'AC240':
        bad_channels = np.hstack((np.arange(0,63), np.arange(number_of_channel-63,number_of_channel)))
    else:
        ValueError('Spectrometer unknown !')
    
    return bad_channels

def compare_Tb_mopi5(self, ds_dict, calibration_cycle):
    fig, axs = plt.subplots(1,1,sharex=True)
    for s in self.spectrometers:
        mask = ds_dict[s].good_channels[calibration_cycle].data
        mask[mask==0]=np.nan
        axs.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,ds_dict[s].Tb[calibration_cycle]*mask, lw=0.5, label=s)
        #axs.set_xlim(110.25, 111.4)
        #axs.set_ylim(0,250)
        #ax].set_ylim(np.median(ds_dict[s].Tb[calibration_cycle].data)-10,np.median(ds_dict[s].Tb[calibration_cycle].data)+15)
        axs.set_xlabel("f [GHz]")
        axs.set_ylabel(r"$T_B$ [K]")
        axs.set_title("Tb")
        axs.grid()
        axs.legend(fontsize='xx-small')
        #ax3.legend()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    return fig

def compare_spectra_mopi5(self, ds_dict, calibration_cycle=0):
    fig, axs = plt.subplots(2,2,sharex=True)
    for s in self.spectrometers:
        mask = ds_dict[s].good_channels[calibration_cycle].data
        mask[mask==0]=np.nan
        axs[0][0].plot(ds_dict[s].frequencies.data/1e9,
                 ds_dict[s].Tb[calibration_cycle].data, lw=0.5, label=s)
        axs[0][0].set_xlim(110.25, 111.4)
        axs[0][0].set_ylim(np.median(ds_dict[s].Tb[calibration_cycle].data)-10,np.median(ds_dict[s].Tb[calibration_cycle].data)+15)
        axs[0][0].set_xlabel("f [GHz]")
        axs[0][0].set_ylabel(r"$T_B$ [K]")
        axs[0][0].set_title("Tb")
        axs[0][0].grid()
        axs[0][0].legend(fontsize='xx-small')
        axs[1][0].plot(ds_dict[s].frequencies.data/1e9,
                 ds_dict[s].Tb_corr_old[calibration_cycle].data*mask,
                 lw=0.5, label=s)
        #axs[1][0].set_xlim(110.3, 111.4)
        axs[1][0].set_ylim(0,30)
        axs[1][0].set_xlabel("f [GHz]")
        axs[1][0].set_ylabel(r"$T_B$ [K]")
        axs[1][0].set_title("Tb_corr_old")
        axs[1][0].grid()
        axs[0][1].plot(ds_dict[s].frequencies.data/1e9,
                 ds_dict[s].Tb_corr[calibration_cycle].data*mask,
                 lw=0.5, label=s)
        #axs[0][1].set_xlim(110.3, 111.4)
        axs[0][1].set_ylim(0,30)
        axs[0][1].set_xlabel("f [GHz]")
        axs[0][1].set_ylabel(r"$T_B$ [K]")
        axs[0][1].set_title("Tb_corr")
        axs[0][1].grid()
        axs[1][1].plot(ds_dict[s].frequencies.data/1e9,
                 ds_dict[s].stdTb[calibration_cycle].data, lw=0.5, label=s)
        #axs[1][1].set_xlim(110.3, 111.4)
        axs[1][1].set_ylim(0,np.median(ds_dict[s].stdTb[calibration_cycle].data)+0.5)
        axs[1][1].set_xlabel("f [GHz]")
        axs[1][1].set_ylabel(r"$stdTb$ [K]")
        axs[1][1].set_title("stdTb")
        axs[1][1].grid()
        #ax3.legend()
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
    
    return fig