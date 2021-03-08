#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Integration script for IAP instruments

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
from abc import ABC
import os
import datetime

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from utils_GROSOM import save_single_pdf

from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.lines import Line2D
import plotly_express as px
from dotenv import load_dotenv

import GROSOM_library



# %%
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.4')
# ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
# ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
# ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

# if __name__ == "__main__":

instrument_name = "GROMOS"

# date = pd.date_range(start='2019-01-03', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

# date = pd.date_range(start='2019-01-30', end='2019-06-18')

date = pd.date_range(start='2019-02-27', end='2019-02-27')
#date = datetime.date(2019,1,15)
#date = [datetime.date(2019,3,11), datetime.date(2019,4,3)]

int_time = 1

plot_cycle = False
df_bins = 200e3

plot_all = False
plot_o3_ts = True
plot_selected = True

integration_strategy = 'classic'
classic = np.arange(1, 24)

cycle = 14
# %%
if instrument_name=="GROMOS":
    import gromos_classes as gc
    basename_lvl1 = "/home/eric/Documents/PhD/GROSOM/Data/"
    basename_lvl2 = "/home/eric/Documents/PhD/GROSOM/Data/"  
    instrument = gc.GROMOS_LvL2(
        date = date, 
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name=="SOMORA":
    import somora_classes as sm
    instrument = sm.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )


 # Plotting part
level2_dataset = instrument.read_level2(
    spectrometers = ['AC240'],
    extra_base='_23'
    )
F0 =  instrument.observation_frequency

if plot_all:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2+instrument.datestr +'_plot_all_test'
    GROSOM_library.plot_O3_all(level2_dataset, outname)

if plot_o3_ts:
    ozone = level2_dataset['AC240'].isel(o3_lat=0, o3_lon=0)
    o3 = ozone.o3_x
    o3 = o3.assign_coords({'altitude':ozone.o3_z})
    #o3=o3.swap_dims({'o3_p':'altitude'})
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
    o3.plot(x='time', y='altitude')


if plot_selected:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2+instrument.datestr +'_plot_sel_test'
    GROSOM_library.plot_O3_all(
        level2_dataset,
        outname,
        cycles=[1,6,12,18]
    )

if plot_cycle:
    figures = []
    level2_dataset = level2_dataset['AC240']

    title = 'Ozone retrievals'
    level2_cycle = level2_dataset.isel(observation=0)

    ozone_ds = level2_cycle.isel(o3_lat=0, o3_lon=0)
    h2o_ds = level2_cycle.h2o.isel(h2o_lat=0, h2o_lon=0)

    print('fshift fit: {:g} kHz'.format(level2_cycle.freq_shift_x[0].values/1e3))
    
    f_backend = level2_cycle.f.values
    y = level2_cycle.y.values
    y_baseline = level2_cycle.y_baseline.values
    yf = level2_cycle.yf.values
    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")
    
    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - F0) / 1e6, y, label="observed")
    axs[0].plot((f_backend - F0) / 1e6, yf, label="fitted")
    #axs[0].set_ylim(-5, 50)
    axs[0].legend()
    if level2_cycle.freq_shift_x[0].values is not None:
        axs[0].text(
            0.02,
            0.8,
            'fshift fit: {:g} kHz'.format(level2_cycle.freq_shift_x[0].values/1e3),
            transform=axs[0].transAxes,
            verticalalignment="bottom",
            horizontalalignment="left",
    )

    axs[1].plot((f_backend - F0) / 1e6, r, label="residuals")
    axs[1].plot((f_backend - F0) / 1e6, r_smooth, label="residuals smooth")
    axs[1].plot((f_backend - F0) / 1e6, y_baseline, label="baseline")
    #axs[1].set_ylim(-2, 2)
    axs[1].legend()
    axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(F0 / 1e9))

    for ax in axs:
        ax.set_ylabel("$T_B$ [K]")
        ax.set_xlim([min((f_backend - F0) / 1e6), max((f_backend - F0) / 1e6)])
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    fig, axs = plt.subplots(1, 2, sharey=True)
    axs[0].plot(
        ozone_ds.o3_x * 1e6, ozone_ds.o3_z / 1e3, label="retrieved", marker="x"
    )
    axs[0].plot(ozone_ds.o3_xa * 1e6, ozone_ds.o3_z / 1e3, label="apriori")
    axs[0].set_xlim(-2,9)
    axs[0].set_xlabel("Ozone VMR [ppm]")
    axs[0].set_ylabel("Altitude [km]")
    axs[0].legend()
    axs[1].plot(ozone_ds.o3_mr, ozone_ds.o3_z / 1e3)
    axs[1].set_xlabel("Measurement response")

    axs[0].grid(True)
    axs[1].grid(True)
    figures.append(fig)

    fig, axs = plt.subplots(1, 2, sharey=True)    

    axs[0].plot(ozone_ds.o3_es * 1e6, ozone_ds.o3_z / 1e3, label="smoothing error")
    axs[0].plot(ozone_ds.o3_eo * 1e6, ozone_ds.o3_z / 1e3, label="obs error")
    axs[0].set_xlabel("$e$ [ppm]")
    axs[0].set_ylabel("Altitude [km]")
    axs[0].legend()

    # axs[1].plot(100*(ozone_ret.x - og_ozone)/og_ozone, ozone_ret.z_grid / 1e3, label="retrieval-og")
    # axs[1].plot(100*(ozone_ret.xa - og_ozone)/og_ozone, z_og / 1e3, label="apriori-og")
    # axs[0].set_xlabel("Rel diff [%]")
    # axs[0].set_ylabel("Altitude [km]")

    for avk in ozone_ds.o3_avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1].plot(avk, ozone_ds.o3_z / 1e3)
    
    axs[1].set_xlabel("AVKM")
    axs[1].grid(True)
    axs[1].grid(True)

    fig.suptitle(title + " Ozone")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)




# # %%
# if plot_bias:
#     end_dates = [datetime.date(2019, 1, 5), datetime.date(
#         2019, 2, 22), datetime.date(2019, 4, 27)]
#     monthly_color = ['magenta', 'blue', 'cyan', 'orange', 'red']
#     month_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun']
#     size = 12
#     figures = list()
#     fig = plt.figure(figsize=(9, 6))
#     ax1 = fig.add_subplot(1, 3, 1)
#     ax2 = fig.add_subplot(1, 3, 2)
#     ax3 = fig.add_subplot(1, 3, 3)
#     # fig2.subplots_adjust(hspace=0.4, wspace=0.1)
#     count = 0
#     symbol = {'U5303': '*', 'AC240': 's'}
#     for d in end_dates:
#         try:
#             integration = mc.MOPI5_LvL2(
#                 d, basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)
#             integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
#                 no_flag=True, meta_data=False)
#             color = monthly_color[d.month-1]
#             for s in ['U5303', 'AC240']:
#                 scatter = ax1.scatter(
#                     integrated_data[s].mean_Tb.data, integrated_data[s].line_amplitude.data, s=size, marker=symbol[s], color=monthly_color[d.month-1])
#                 ax1.set_ylabel(r'$T_B$ [K]')
#                 # ax1.set_ylim(-2,0)
#                 # ax1.set_ylim(-1,0.5)
#                 ax1.set_title('line amplitude')
#                 ax1.set_xlabel('Mean $T_b$ [K]')
#                 ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].continuum_value_line_center.data,
#                             s=size, marker=symbol[s], color=monthly_color[d.month-1])
#                 ax2.set_ylabel(r'$T_B$ [K]')
#                 ax2.set_xlabel('Mean $T_b$ [K]')
#                 ax2.set_title('$T_B$ continuum')
#                 # ax2.set_ylim(-2,0)
#                 # ax2.set_ylim(-1,1)
#                 ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_indiv *
#                             1e9, s=size, marker=symbol[s], color=monthly_color[d.month-1])
#                 ax3.set_ylabel('m [K/GHz]')
#                 # ax3.set_ylim(-0.8,0.2)
#                 ax3.set_title('Slope')
#                 ax3.set_xlabel('Mean $T_b$ [K]')
         
#         except: 
#             print('no data for :', d)
#             pass
#     legend_elements_both = [
#         Line2D([0], [0], marker=symbol['U5303'], color='w', markerfacecolor='magenta',
#                label='U5303 : '+month_name[0], markersize=size-2),
#         Line2D([0], [0], marker=symbol['U5303'], color='w', markerfacecolor='blue',
#                label='U5303 : '+month_name[1], markersize=size-2),
#         Line2D([0], [0], marker=symbol['U5303'], color='w', markerfacecolor='orange',
#                label='U5303 : '+month_name[3], markersize=size-2),
#         Line2D([0], [0], marker=symbol['AC240'], color='w', markerfacecolor='magenta',
#                label='AC240 : '+month_name[0], markersize=size-2),
#         Line2D([0], [0], marker=symbol['AC240'], color='w', markerfacecolor='blue',
#                label='AC240 : '+month_name[1], markersize=size-2),
#         Line2D([0], [0], marker=symbol['AC240'], color='w', markerfacecolor='orange',
#                label='AC240 : '+month_name[3], markersize=size-2)
#     ]
#     # legend = ax2.legend(*scatter.legend_elements(prop='marker'),['U5303','AC240'] , fontsize='xx-small',loc=1, title='Month')
#     # ax2.add_artist(legend)
#     # ax2.legend(['U5303','AC240'], fontsize='small')
#     ax1.legend(handles=legend_elements_both, fontsize='small', loc=1)
#     # ax2.legend(handles=legend_elements_spectro, fontsize='small',loc=1)
#     ax1.grid()
#     ax2.grid()
#     ax3.grid()
#     fig.suptitle('Absolute bias')
#     fig.tight_layout(rect=[0, 0.01, 1, 0.95])
#     plt.show()

# %%
