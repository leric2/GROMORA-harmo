#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib.lines import Line2D
from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator)

from abc import ABC
import os
import datetime

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from utils_GROSOM import save_single_pdf
from dotenv import load_dotenv
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]})

plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 17

load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.4')

# cmpa_str = 'berlin'  # batlow, devon, oslo, imola, lapaz
# cm_data = np.loadtxt(
#     '/home/eric/Documents/PhD/ScientificColourMaps7/'+cmpa_str+'/'+cmpa_str+'.txt')
# cmap_crameri = LinearSegmentedColormap.from_list('berlin', cm_data)

import mopi5_library
import mopi5_classes as mc
# %%

# if __name__ == "__main__":

instrument_name = "mopi5"
# date = datetime.date(2019,2,21)
# date = pd.date_range(start='2019-01-05', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

# 
# date = pd.date_range(start='2019-01-30', end='2019-06-18')

date = pd.date_range(start='2019-02-22', end='2019-02-22')
meanTb_chunks = [80, 85, 90, 95, 100, 105,
                 110, 115, 120, 130, 140, 150, 170, 190]
lowerBound = [0, 80, 85, 90, 95, 100, 105,
              110, 115, 120, 130, 140, 150, 170, 190]

# date = pd.date_range(start='2019-05-01', end='2019-05-04')
# No U5303

# date = pd.date_range(start='2019-04-27', end='2019-04-27')
# meanTb_chunks = [105, 110, 115, 120, 130, 160, 180, 200]
# lowerBound = [0, 105, 110, 115, 120, 130, 160, 180, 200]

# date = pd.date_range(start='2019-06-15', end='2019-06-15')
# meanTb_chunks = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220]
# lowerBound = [0, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220]

# date = pd.date_range(start='2019-03-12', end='2019-03-12')
# date = pd.date_range(start='2019-02-22', end='2019-02-22')
# options are: 'TOD', 'TOD_harmo', 'classic' 'meanTb_harmo', or 'meanTb'
integration_strategy = 'meanTb_harmo'
int_time = 1

plot_ts_Tb_Tsys = False
df_bins = 200e3
#date1b = pd.to_datetime(date[-1])

plot_spectra_schematic = False

plot_comparison = False
compare_level2_mopi5 = True
compare_alpha = False

plot_spectra_comparison_scaling_corr_paper = False
plot_spectra_comparison_3_spectro_paper = True
plot_bias = False
plot_bias_TOD = False
plot_bias_TOD_full = False
plot_o3 = False
plot_o3_sel = False
plot_sel_paper = False
plot_avks_paper = False
plot_non_lin = False

plot_bias_spectra_monthly = False

# Define the parameters for integration
# TOD = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]
# interval = np.ones(len(TOD))
TOD = np.arange(24)
interval = 0.5*np.ones(len(TOD))
# TOD = [3, 9, 15, 21]

classic = np.arange(1, 24)

# %%
outfolder = '/home/eric/Documents/PhD/MOPI/Data/Level3/'
#outfolder = '/home/esauvageat/Documents/MOPI5/Level3/'
basename_lvl1 = "/storage/tub/instruments/mopi5/level1/"
basename_lvl2 = "/scratch/MOPI5/Level1/"
basename_lvl2 = "/storage/tub/instruments/mopi5/level2/"
#basename_lvl1 = "/home/eric/Documents/PhD/MOPI/Data/Level1b_test_Gerald/"
#basename_lvl2 = "/home/eric/Documents/PhD/MOPI/Data/Level2/"
# calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time, ['AC240','USRP-A'])
integration = mc.MOPI5_LvL2(date[0], basename_lvl1, basename_lvl2,
                            integration_strategy, integration_time=int_time)
# Plotting part
integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
    no_flag=True, meta_data=False)
if integration_strategy == 'meanTb' or integration_strategy == 'meanTb_harmo':
    identifier_plot = integrated_data['AC240'].coords['chunks'].data[1:].tolist(
    ) + [300]
    # identifier_plot = meanTb_chunks  + [300]
    idx_all = np.arange(0, len(identifier_plot))
    dimension = ['chunks', 'channel_idx']
else:
    dimension = ['time', 'channel_idx']
    identifier_plot = TOD
    idx_all = np.arange(0, len(TOD))
param_slope_broadband = [111.2e9, 25e6, 110.5e9, 25e6]
param_slope = {'AC240': param_slope_broadband,
               'USRP-A': [110.84e9, 5e6, 110.72e9, 5e6], 'U5303': param_slope_broadband}

theoretical_nonlinearities = np.polyfit(
        [80, 186, 292], [0, -0.20, 0], deg=2)
fitted_poly_theoretical_nonlinearities = np.poly1d(
        theoretical_nonlinearities)

if plot_spectra_schematic:
    integration.compare_spectra_mopi5(
        spectrometers=['AC240'],
        dim=dimension[0],
        idx=idx_all,
        # idx=[0,1,2,3],
        save_plot=True,
        identifier=identifier_plot,
        lowerBound=lowerBound,
        with_corr=False,
        corr_band=param_slope,
        title='',
    )
    
# %%

if plot_spectra_comparison_3_spectro_paper:
    theoretical_nonlinearities = np.polyfit(
        [80, 186, 292], [0, -0.20, 0], deg=2)
    fitted_poly_theoretical_nonlinearities = np.poly1d(theoretical_nonlinearities)
    f_bands = 1e6*np.arange(-450,550,100)
    f_band_width = 100e6

    # Mean comparison
    df_around_line = 10e6
    mean_bias_USRP = np.ones(len(idx_all))
    for i in idx_all:
        Tb_diff = integrated_data['USRP-A'].interpolated_Tb[i].data - integrated_data['U5303'].binned_Tb[i].data
        Tb_diff = Tb_diff[(integrated_data['USRP-A'].bin_freq_interp> integration.observation_frequency-df_around_line ) & (integrated_data['USRP-A'].bin_freq_interp< integration.observation_frequency+df_around_line )]
        # print('Mean bias around obs freq for USRP-A an cycle ',str(i),': ')
        # print(np.nanmean(Tb_diff))
        mean_bias_USRP[i] = np.nanmean(Tb_diff)
    
    
    mean_bias_AC240 = np.ones(len(idx_all))
    mean_Tb = np.ones(len(idx_all))
    for i in idx_all:
        Tb_diff = integrated_data['AC240'].interpolated_Tb[i].data - integrated_data['U5303'].binned_Tb[i].data
        Tb_diff = Tb_diff[(integrated_data['AC240'].bin_freq_interp> integration.observation_frequency-df_around_line ) & (integrated_data['AC240'].bin_freq_interp< integration.observation_frequency+df_around_line )]
        # print('Mean bias around obs freq for AC240 an cycle ',str(i),': ')
        # print(np.nanmean(Tb_diff))
        mean_bias_AC240[i] = np.nanmean(Tb_diff)
        mean_Tb[i]= np.nanmean(integrated_data['U5303'].binned_Tb[i].data)

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(mean_Tb, mean_bias_AC240, '-o')
    ax.grid()
    ax.set_title(r'Differences AC240-U5303 around line center ($\pm$ 10 MHz)')
    ax.set_xlabel('Tb [K]',fontsize=12)
    ax.set_ylabel('Tb differences [K]',fontsize=12)
    fig.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig.savefig('/storage/nas/MW/scratch/sauvageat/MOPI5/ToDiscuss/' + 'Tb_diff_line_center.pdf')

    rel_bias_freq_bands =  np.ones((len(idx_all), len(f_bands)))
    mean_tb_bands = np.ones((len(idx_all), len(f_bands)))
    a=0.08
    for i in idx_all:
        Tb_corr = (1/(1-a))*(integrated_data['AC240'].interpolated_Tb[i].data-a*np.nanmean(integrated_data['AC240'].interpolated_Tb[i].data))
        #Tb_diff = integrated_data['AC240'].interpolated_Tb[i].data - integrated_data['U5303'].binned_Tb[i].data
        Tb_diff = Tb_corr - integrated_data['U5303'].binned_Tb[i].data
        Tb_ref = integrated_data['U5303'].binned_Tb[i].data
        for j, f in enumerate(f_bands):
            f_start = integration.observation_frequency+f-f_band_width/2
            f_stop = integration.observation_frequency+f+f_band_width/2
            diff = Tb_diff[(integrated_data['AC240'].bin_freq_interp> f_start ) & (integrated_data['AC240'].bin_freq_interp< f_stop )]
            Tb_mean_band = Tb_ref[(integrated_data['AC240'].bin_freq_interp> f_start ) & (integrated_data['AC240'].bin_freq_interp< f_stop )]
            rel_bias_freq_bands[i, j] = np.nanmean(diff)#/np.nanmean(Tb_mean_band)
            mean_tb_bands[i, j] = np.nanmean(Tb_mean_band)
    
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    for m in range(len(idx_all)):
        ax.plot(1e-9*(integration.observation_frequency+f_bands), rel_bias_freq_bands[m,:], '-o', label='Tb = '+str(identifier_plot[m])+' K')
    ax.grid()
    ax.legend(fontsize='small',loc='upper right', ncol=2)
    ax.set_title('Differences AC240-U5303')
    ax.set_xlabel('Frequency [GHz]',fontsize=12)
    ax.set_ylabel('Tb differences [K]',fontsize=12)
    fig.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig.savefig('/storage/nas/MW/scratch/sauvageat/MOPI5/ToDiscuss/' + 'Tb_diff_vs_freq_bands.pdf')

    fig, axes = plt.subplots(nrows=2, ncols=5,figsize=(20, 15))
    for fn in [0,1,2,3,4]:
        axes[0,fn].plot(mean_tb_bands[:,fn], rel_bias_freq_bands[:,fn], '-o', label='Non-linearity term')
        axes[0,fn].plot(np.arange(75,220,2), fitted_poly_theoretical_nonlinearities(np.arange(75,220,2)), label='broadband correction')
        axes[0,fn].grid()
        axes[0,fn].set_ylim(-0.35, 0.1)
        axes[0,fn].set_title(f'f = {1e-9*(integration.observation_frequency+f_bands[fn]):.2f} GHz' )
        fitted_non_lin_per_frequency_band = np.polyfit(mean_tb_bands[:,fn], rel_bias_freq_bands[:,fn], deg=2)
        print(fitted_non_lin_per_frequency_band[0])
        fitted_non_lin_freq_band_plot = np.poly1d(fitted_non_lin_per_frequency_band)
        axes[0,fn].plot(np.arange(75,220,2), fitted_non_lin_freq_band_plot(np.arange(75,220,2)), label=f'fit, QuadTerm = {1e5*fitted_non_lin_per_frequency_band[0]:.1f}e-5')
        axes[0,fn].legend()
    for fn in [5,6,7,8,9]:
        axes[1,fn-5].plot(mean_tb_bands[:,fn], rel_bias_freq_bands[:,fn], '-o', label='Non-linearity term')
        axes[1,fn-5].plot(np.arange(75,220,2), fitted_poly_theoretical_nonlinearities(np.arange(75,220,2)), label='broadband correction')
        axes[1,fn-5].grid()
        axes[1,fn-5].set_xlabel('Tb [K]',fontsize=12)
        axes[1,fn-5].set_ylim(-0.35, 0.1)
        axes[1,fn-5].set_title(f'f = {1e-9*(integration.observation_frequency+f_bands[fn]):.2f} GHz' )
        fitted_non_lin_per_frequency_band = np.polyfit(mean_tb_bands[:,fn], rel_bias_freq_bands[:,fn], deg=2)
        print(fitted_non_lin_per_frequency_band[0])
        fitted_non_lin_freq_band_plot = np.poly1d(fitted_non_lin_per_frequency_band)
        axes[1,fn-5].plot(np.arange(75,220,2), fitted_non_lin_freq_band_plot(np.arange(75,220,2)), label=f'fit, QuadTerm = {1e5*fitted_non_lin_per_frequency_band[0]:.1f}e-5')
        axes[1,fn-5].legend()
    #ax.legend(fontsize='small',loc='upper right')
    fig.suptitle('Differences AC240 (unscaled)-U5303 by 100 MHz frequency bands', fontsize=20)
    axes[0,0].set_ylabel('Tb differences [K]',fontsize=16)
    axes[1,0].set_ylabel('Tb differences [K]',fontsize=16)
    fig.tight_layout(rect=[0, 0.01, 1, 0.95])

    fig.savefig('/storage/nas/MW/scratch/sauvageat/MOPI5/ToDiscuss/' + 'Tb_diff_vs_Tb_with_non-lin_theory.pdf')

    mean_bias_AC240_corr = np.ones(len(idx_all))
    a = 0.08
    for i in idx_all:
        Tb_ac240 = integrated_data['AC240'].interpolated_Tb[i].data
        non_lin = fitted_poly_theoretical_nonlinearities(integrated_data['AC240'].interpolated_Tb[i].data)
        Tb_ac240_corr = (1/(1-a))*(Tb_ac240 - a*np.nanmean(Tb_ac240) - non_lin)
        Tb_diff = Tb_ac240_corr - integrated_data['U5303'].binned_Tb[i].data
        Tb_diff = Tb_diff[(integrated_data['AC240'].bin_freq_interp> integration.observation_frequency-df_around_line ) & (integrated_data['AC240'].bin_freq_interp< integration.observation_frequency+df_around_line )]
        print('Mean bias around obs freq for AC240 corrected an cycle ',str(i),': ')
        print(np.nanmean(Tb_diff))
        mean_bias_AC240_corr[i] = np.nanmean(Tb_diff)

    np.savetxt(outfolder + 'meanBiasUSRP.txt',1000*mean_bias_USRP,fmt='%.0f')
    np.savetxt(outfolder + 'meanBiasAC240.txt',1000*mean_bias_AC240,fmt='%.0f')
    np.savetxt(outfolder + 'meanBiasAC240Corr.txt',1000*mean_bias_AC240_corr,fmt='%.0f')

    integrated_data['U5303'].frequencies 
    integration.compare_spectra_binned_interp_mopi5(
        dim=dimension[0],
        idx=[0,11],#np.arange(0,15),
        spectrometers=['AC240', 'USRP-A'],
        save_plot=True,
        use_basis='U5303',
        # identifier=TOD,
        identifier=identifier_plot,
        clean=True,
        corrected=False,
        outfolder=outfolder
    )
if plot_non_lin:
    continuum_amplitude_diff = integrated_data['AC240'].continuum_value_line_center - \
        integrated_data['U5303'].continuum_value_line_center
    fitted_factor_continuum = np.polyfit(
        integrated_data['AC240'].mean_Tb, continuum_amplitude_diff, deg=2)
    fitted_poly_continuum = np.poly1d(fitted_factor_continuum)
    plt.plot(fitted_poly_theoretical_nonlinearities(integrated_data['U5303'].interpolated_Tb[10].data))


if plot_bias_spectra_monthly:
#%%

    #     integration.compare_spectra_binned_interp_mopi5(
    #         dim=dimension[0],
    #         idx=idx_all,
    #         spectrometers=['AC240','USRP-A'],
    #         save_plot = True,
    #         use_basis='U5303',
    #         #identifier=TOD,
    #         identifier=identifier_plot,
    #         clean=True
    # )
    # Fitting the correction factors
    #line_amplitude_diff = (integrated_data['AC240'].line_amplitude -integrated_data['U5303'].line_amplitude)/integrated_data['U5303'].line_amplitude
    continuum_amplitude_diff = integrated_data['AC240'].continuum_value_line_center - \
        integrated_data['U5303'].continuum_value_line_center
    fitted_factor_continuum = np.polyfit(
        integrated_data['AC240'].mean_Tb, continuum_amplitude_diff, deg=2)
    fitted_poly_continuum = np.poly1d(fitted_factor_continuum)

    s = 'U5303'
    # plt.plot(integrated_data[s].bin_freq_interp.data, integrated_data[s].interpolated_Tb[11].data, 'k-', linewidth=0.4)
    # plt.plot(integrated_data[s].bin_freq_interp.data, integrated_data[s].interpolated_Tb[11].data+fitted_poly_theoretical_nonlinearities(integrated_data[s].interpolated_Tb[11].data), 'r-', linewidth=0.4)
    non_lin = fitted_poly_theoretical_nonlinearities(
        integrated_data[s].interpolated_Tb.data)
    param_slope = {'AC240': param_slope_broadband,
                   'USRP-A': [110.84e9, 5e6, 110.72e9, 5e6], 'U5303': param_slope_broadband}
    # integrated_data = integration.add_bias_characterization_with_non_linearities(
    #     spectrometers=['U5303'],
    #     dim=dimension[0],
    #     param_slope=param_slope,
    #     around_center_value=1e6
    # )
    fig1 = plt.figure(figsize=(9, 6))
    ax1 = fig1.add_subplot(1, 3, 1)
    ax2 = fig1.add_subplot(1, 3, 2)
    ax3 = fig1.add_subplot(1, 3, 3)
    color = 'b'
    s = 'AC240'
    scatter = ax1.scatter(integrated_data[s].mean_Tb.data,(integrated_data[s].line_amplitude.data -
                          integrated_data['U5303'].line_amplitude.data), s=10, color=color)
    ax1.set_ylabel(r'$\Delta T_B$ [\%]')
    # ax1.set_ylim(-2,0)
    # ax1.set_ylim(-1,0.5)
    ax1.set_title('line amplitude difference')
    ax1.set_xlabel('Mean $T_B$ [K]')
    # ax1.set_ylim(-0.92,-0.25)
    ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].continuum_value_line_center.data -
                integrated_data['U5303'].continuum_value_line_center.data, s=10, color=color)
    ax2.set_ylabel(r'$\Delta T_B$ [K]')
    ax2.set_xlabel('Mean $T_B$ [K]')
    ax2.set_title(r'$\Delta T_B$ continuum')
    #ax2.set_ylim(-0.22, 0.06)
    # ax2.set_ylim(-1,1)
    ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_indiv *
                1e9-integrated_data['U5303'].slope_indiv.data*1e9, s=10, color=color)
    ax3.set_ylabel(r'$\Delta m$ [K/GHz]')
    # ax3.set_ylim(-0.8,0.2)
    #ax3.set_ylim(-0.55, 0.05)
    ax3.set_title('Slope difference')
    ax3.set_xlabel('Mean $T_B$ [K]')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig1.suptitle('Difference AC240 - U5303')
    fig1.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig1.savefig('/home/eric/Documents/PhD/MOPI/Data/Level3/' + 'bias_June.pdf')
    fig2 = plt.figure(figsize=(9, 6))
    ax1 = fig2.add_subplot(1, 3, 1)
    ax2 = fig2.add_subplot(1, 3, 2)
    ax3 = fig2.add_subplot(1, 3, 3)
    color = 'b'
    s = 'AC240'
    scatter = ax1.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].line_amplitude.data -
                          integrated_data['U5303'].line_amplitude.data)/integrated_data['U5303'].line_amplitude.data, s=10, color=color)
    ax1.set_ylabel(r'$\Delta T_B$ [\%]')
    # ax1.set_ylim(-2,0)
    # ax1.set_ylim(-1,0.5)
    ax1.set_title('line amplitude difference')
    ax1.set_xlabel('Mean $T_B$ [K]')
    # ax1.set_ylim(-0.92,-0.25)
    ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].continuum_value_line_center.data -
                integrated_data['U5303'].continuum_value_line_center.data, s=10, color=color)
    ax2.set_ylabel(r'$\Delta T_B$ [K]')
    ax2.set_xlabel('Mean $T_B$ [K]')
    ax2.set_title(r'$\Delta T_B$ continuum')
    #ax2.set_ylim(-0.22, 0.06)
    # ax2.set_ylim(-1,1)
    ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_indiv *
                1e9-integrated_data['U5303'].slope_indiv.data*1e9, s=10, color=color)
    ax3.set_ylabel(r'$\Delta m$ [K/GHz]')
    # ax3.set_ylim(-0.8,0.2)
    #ax3.set_ylim(-0.55, 0.05)
    ax3.set_title('Slope difference')
    ax3.set_xlabel('Mean $T_B$ [K]')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig2.suptitle('Difference AC240 - U5303')
    fig2.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig2.savefig('/home/eric/Documents/PhD/MOPI/Data/Level3/' + 'bias_relative_June.pdf')
    plt.show()

    # plt.plot(integrated_data['AC240'].bin_freq_interp.data/1e9, non_lin[0])
    # plt.plot(integrated_data['AC240'].bin_freq_interp.data/1e9, non_lin[11])
    # plt.ylabel(r'$\Delta T_B$ [K] nonlin')
    # plt.xlabel('Frequency [GHz]')
    # plt.legend((r'$T_B < 80K$', r'$T_B >140K$'))


if plot_spectra_comparison_scaling_corr_paper:
    integration.compare_spectra_binned_interp_mopi5_factor(
        dim=dimension[0],
        idx=[0,11],
        spectrometers=['AC240'],
        save_plot=True,
        use_basis='U5303',
        # identifier=TOD,
        identifier=identifier_plot,
        # alpha=-100*line_amplitude_diff.data,
        alpha=8*np.ones(15),
        binning=2,
        lowerBound=lowerBound,
        variable=True,
        broadband_bias=fitted_poly_theoretical_nonlinearities(integrated_data['AC240'].interpolated_Tb.data),
        paper=True,
        outfolder=outfolder
    )

if compare_alpha:
    integration.compare_spectra_binned_interp_mopi5_factor(
        dim=dimension[0],
        idx=[0,11],
        spectrometers=['AC240'],
        save_plot = True,
        use_basis='U5303',
        #identifier=TOD,
        identifier=identifier_plot,
        #alpha=8*np.ones(15),
        alpha=[7,8,9],
        binning=4,
        lowerBound=lowerBound,
        variable=True,
        broadband_bias=fitted_poly_theoretical_nonlinearities(integrated_data['AC240'].interpolated_Tb.data),
        outfolder=outfolder
    )

    # %%
if plot_comparison:
    # integration.compare_spectra_mopi5(
    #     dim=dimension[0],
    #     idx=idx_all,
    #     # idx=[0,1,2,3],
    #     save_plot=True,
    #     # identifier=TOD,
    #     identifier=identifier_plot,
    #     with_corr=True,
    #     corr_band=False
    # )
    integration.plot_time_min_comp()
    integration.compare_spectra_binned_interp_mopi5(
        dim=dimension[0],
        idx=idx_all,
        # idx=[0,1,2,3],
        spectrometers=integration.spectrometers,
        save_plot=True,
        # identifier=TOD,
        identifier=identifier_plot,
        use_basis='U5303'
    )
    integration.compare_spectra_binned_interp_mopi5(
        dim=dimension[0],
        # idx=np.arange(0,len(meanTb_chunks)+1),
        idx=idx_all,
        spectrometers=integration.spectrometers,
        save_plot=True,
        identifier=identifier_plot,
        use_basis='U5303',
        corrected=True
    )

# %%
if plot_bias:
    # end_dates = [datetime.date(2019, 1, 5), datetime.date(
    #     2019, 2, 22), datetime.date(2019, 4, 27)]
    end_dates = [pd.date_range(start='2019-02-22', end='2019-02-22')]
    monthly_color = ['magenta', 'blue', 'cyan', 'orange', 'red']
    month_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun']
    size = 12
    figures = list()
    fig = plt.figure(figsize=(9, 6))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)
    # fig2.subplots_adjust(hspace=0.4, wspace=0.1)
    count = 0
    symbol = {'U5303': '*', 'AC240': 's'}
    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(
                d, basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
                no_flag=True, meta_data=False)
            #integrated_data = integration.add_mean_Tb(spectrometers = integration.spectrometers, around_center=True, around_center_value=500e6)

            color = monthly_color[d[0].month-1]
            for s in ['U5303', 'AC240']:
                scatter = ax1.scatter(
                    integrated_data[s].mean_Tb.data, integrated_data[s].line_amplitude.data, s=size, marker=symbol[s], color=monthly_color[d[0].month-1])
                ax1.set_ylabel(r'$T_B$ [K]')
                # ax1.set_ylim(-2,0)
                # ax1.set_ylim(-1,0.5)
                ax1.set_title('line amplitude')
                ax1.set_xlabel('Mean $T_b$ [K]')
                ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].continuum_value_line_center.data,
                            s=size, marker=symbol[s], color=monthly_color[d[0].month-1])
                ax2.set_ylabel(r'$T_B$ [K]')
                ax2.set_xlabel('Mean $T_b$ [K]')
                ax2.set_title('$T_B$ continuum')
                # ax2.set_ylim(-2,0)
                # ax2.set_ylim(-1,1)
                ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_indiv *
                            1e9, s=size, marker=symbol[s], color=monthly_color[d[0].month-1])
                ax3.set_ylabel('m [K/GHz]')
                # ax3.set_ylim(-0.8,0.2)
                ax3.set_title('Slope')
                ax3.set_xlabel('Mean $T_b$ [K]')

        except:
            print('no data for :', d)
            pass
    legend_elements_both = [
        Line2D([0], [0], marker=symbol['U5303'], color='w', markerfacecolor='magenta',
               label='U5303 : '+month_name[0], markersize=size-2),
        Line2D([0], [0], marker=symbol['U5303'], color='w', markerfacecolor='blue',
               label='U5303 : '+month_name[1], markersize=size-2),
        Line2D([0], [0], marker=symbol['U5303'], color='w', markerfacecolor='orange',
               label='U5303 : '+month_name[3], markersize=size-2),
        Line2D([0], [0], marker=symbol['AC240'], color='w', markerfacecolor='magenta',
               label='AC240 : '+month_name[0], markersize=size-2),
        Line2D([0], [0], marker=symbol['AC240'], color='w', markerfacecolor='blue',
               label='AC240 : '+month_name[1], markersize=size-2),
        Line2D([0], [0], marker=symbol['AC240'], color='w', markerfacecolor='orange',
               label='AC240 : '+month_name[3], markersize=size-2)
    ]
    # legend = ax2.legend(*scatter.legend_elements(prop='marker'),['U5303','AC240'] , fontsize='xx-small',loc=1, title='Month')
    # ax2.add_artist(legend)
    # ax2.legend(['U5303','AC240'], fontsize='small')
    ax1.legend(handles=legend_elements_both, fontsize='small', loc=1)
    # ax2.legend(handles=legend_elements_spectro, fontsize='small',loc=1)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig.suptitle('Absolute bias')
    fig.tight_layout(rect=[0, 0.01, 1, 0.95])
    plt.show()
    figures.append(fig)
    # save_single_pdf(basename_lvl1+'full_bias_-500-200Mhz'+'_'+integration_strategy+'.pdf',figures)
    fig2 = plt.figure(figsize=(9, 6))
    ax1 = fig2.add_subplot(1, 3, 1)
    ax2 = fig2.add_subplot(1, 3, 2)
    ax3 = fig2.add_subplot(1, 3, 3)
    fig2.subplots_adjust(hspace=0.1, wspace=0)
    count = 0
    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(
                d, basename_lvl1, integration_strategy, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
                no_flag=True, meta_data=False)
            color = monthly_color[d[0].month-1]
            s = 'AC240'
            scatter = ax1.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].line_amplitude.data -
                                  integrated_data['U5303'].line_amplitude.data, s=size, color=color)
            ax1.set_ylabel(r'$\Delta T_B$ [K]')
            # ax1.set_ylim(-2,0)
            # ax1.set_ylim(-1,0.5)
            ax1.set_title('line amplitude difference')
            ax1.set_xlabel('Mean $T_B$ [K]')
            THot = np.mean(integrated_data[s].THot.data)
            theoretical_nonlinearities = np.polyfit(
                [80, (THot+80)/2, THot], [0, -0.2, 0], deg=2)
            fitted_poly_theoretical_nonlinearities = np.poly1d(
                theoretical_nonlinearities)
            ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].continuum_value_line_center.data -
                        integrated_data['U5303'].continuum_value_line_center.data, s=size, color=color)
            ax2.plot(np.arange(70, 310, 10), fitted_poly_theoretical_nonlinearities(
                np.arange(70, 310, 10)), '-', color=color)
            ax2.set_ylabel(r'$\Delta T_B$ [K]')
            ax2.set_xlabel('Mean $T_B$ [K]')
            ax2.set_title(r'$\Delta T_B$ continuum')
            # ax2.set_ylim(-2,0)
            # ax2.set_ylim(-1,1)
            ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_indiv *
                        1e9-integrated_data['U5303'].slope_indiv.data*1e9, s=size, color=color)
            ax3.set_ylabel(r'$\Delta m$ [K/GHz]')
            # ax3.set_ylim(-0.8,0.2)
            ax3.set_title('Slope difference')
            ax3.set_xlabel('Mean $T_B$ [K]')
        except:
            print('no data for :', d)
            pass
    legend_elements = [
        Line2D([0], [0], marker='.', color='w', markerfacecolor='magenta',
               label=month_name[0], markersize=size+2),
        Line2D([0], [0], marker='.', color='w', markerfacecolor='blue',
               label=month_name[1], markersize=size+2),
        Line2D([0], [0], marker='.', color='w', markerfacecolor='orange',
               label=month_name[3], markersize=size+2)
    ]
    # legend = ax1.legend(*scatter.legend_elements(prop='colors'), month_name, fontsize='xx-small',loc=1, title='Month')
    # ax1.add_artist(legend)
    # ax2.legend(['U5303','AC240'], fontsize='small')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig2.suptitle('Difference : '+s+' - U5303')
    ax1.legend(handles=legend_elements, fontsize='small', loc=1)
    # fig.suptitle('Mean hot counts')
    fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    figures.append(fig2)
    fig3 = plt.figure(figsize=(6, 9))
    ax3 = fig3.add_subplot(1, 3, 3)
    ax1 = fig3.add_subplot(1, 3, 1)
    ax2 = fig3.add_subplot(1, 3, 2)
    size = 6
    count = 0
    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(
                d, basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
                no_flag=True, meta_data=False)
            color = monthly_color[d[0].month-1]
            s = 'AC240'
            scatter = ax1.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].line_amplitude.data -
                                                                        integrated_data['U5303'].line_amplitude.data)/integrated_data['U5303'].line_amplitude.data, color=color, s=12)
            # scatter = ax1.scatter(integrated_data[s].time_of_day.data, integrated_data[s].THot.data-integrated_data['U5303'].line_amplitude.data, color=color, s=12)
            ax1.set_ylabel(r'$\Delta T_B$ [%]')
            # ax1.set_ylim(-2,0)
            # ax1.set_ylim(-1,0.5)
            ax1.set_title('line amplitude difference')
            ax1.set_xlabel('Mean $T_B$ [K]')
            ax2.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].continuum_value_line_center.data -
                                                              integrated_data['U5303'].continuum_value_line_center.data)/integrated_data['U5303'].continuum_value_line_center.data, color=color, s=12)
            ax2.set_ylabel(r'$\Delta T_B$ [%]')
            ax2.set_xlabel('Mean $T_B$ [K]')
            ax2.set_title(r'$\Delta T_B$ continuum')
            # ax2.set_ylim(-2,0)
            # ax2.set_ylim(-1,1)
            ax3.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].slope_indiv -
                                                              integrated_data['U5303'].slope_indiv)/integrated_data['U5303'].slope_indiv, color=color, s=12)
            ax3.set_ylabel(r'$\Delta m$ [%]')
            # ax3.set_ylim(-0.8,0.2)
            ax3.set_title('Slope difference')
            ax3.set_xlabel('Mean $T_B$ [K]')
        except:
            print('no data for :', d)
            pass
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax3.legend(handles=legend_elements, fontsize='small', loc=1)
    fig3.suptitle('Fractional Difference : '+s+' - U5303')
    # fig.suptitle('Mean hot counts')
    fig3.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    figures.append(fig3)
    # save_single_pdf(basename_lvl1+'full_bias_feb'+s+'_'+integration_strategy+'.pdf',figures)
    save_single_pdf(basename_lvl1+'full_bias_diff_all' +
                    '_'+integration_strategy+'.pdf', figures)
    figures3 = list()
    fig4 = plt.figure(figsize=(5, 5))
    ax1 = fig4.add_subplot(1, 1, 1)
    # ax2 = fig4.add_subplot(1,3,2)
    # ax3 = fig4.add_subplot(1,3,3)
    count = 0
    size = 5
    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(
                d, basename_lvl1, integration_strategy, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
                no_flag=True, meta_data=False)
            color = monthly_color[d[0].month-1]
            ax1.scatter(integrated_data['AC240'].mean_Tb.data, 100 *
                        integrated_data['AC240'].bias_Tb_lc_corr.data, marker='x', color=color, s=10)
            ax1.scatter(integrated_data['USRP-A'].mean_Tb.data, 100 *
                        integrated_data['USRP-A'].bias_Tb_lc_corr.data, marker='o', color=color, s=10)
            ax1.set_ylabel(r'$\Delta T_B$ [%]')
            # ax1.set_ylim(-2,0)
            # ax1.set_ylim(-1,0.5)
            ax1.set_title('Line center bias after tropospheric correction')
            ax1.set_xlabel('Mean $T_B$ [K]')
            # ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_corr.data*1e9, color=color, s=size)
            # ax2.set_ylabel(r'$\Delta T_B$ [K]')
            # ax2.set_xlabel('Mean $T_B$ [K]')
            # ax2.set_title('Corrected slope continuum')
            # ax2.set_ylim(-2,0)
            # ax2.set_ylim(-1,1)
            # ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_indiv*1e9-integrated_data['U5303'].slope_indiv.data*1e9, color=color, s=size)
            # ax3.set_ylabel('$\Delta m$ [K/GHz]')
            # #ax3.set_ylim(-0.8,0.2)
            # ax3.set_title('Slope difference')
            # ax3.set_xlabel('Mean $T_B$ [K]')
        except:
            print('no data for :', d)
            pass
    # legend = ax1.legend(*scatter.legend_elements(prop='colors'), month_name, fontsize='xx-small',loc=1, title='Month')
    # ax1.add_artist(legend)
    # ax2.legend(['U5303','AC240'], fontsize='small')
    ax1.text(80, -8, 'AC240', fontsize=12)
    ax1.text(180, -1.1, 'USRP-A', fontsize=12)
    ax1.grid()
    # ax2.grid()
    ax1.legend(handles=legend_elements, fontsize=10, loc=5, title='Month')
    fig4.suptitle('Fractional difference with U5303')
    # fig.suptitle('Mean hot counts')
    fig4.tight_layout(rect=[0, 0.01, 1, 0.95])
    plt.show()
    figures3.append(fig4)
    save_single_pdf(basename_lvl1+'bias_lc_tropoCorr' +
                    '_'+integration_strategy+'.pdf', figures3)

# %%
if plot_bias_TOD:
    import matplotlib
    cmap = matplotlib.cm.get_cmap('Dark2') 
    #monthly_color = ['indigo', 'green', 'darkorange', 'gold', 'red']
    monthly_color = [cmap(0.01), cmap(0.28),cmap(0.97), cmap(0.63),  'red']
    integration_strategy = 'TOD_harmo'
    month_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun']
    symbols = ['v','<','^','>','v']
    size = 8
    figures2 = list()
    end_dates = pd.date_range(start='2019-01-03', end='2019-04-30')
    fig5a = plt.figure(figsize=(4,4))
    fig5b = plt.figure(figsize=(4,4))
    fig5c = plt.figure(figsize=(4,4))
    ax1 = fig5a.add_subplot(1, 1, 1)
    ax2 = fig5b.add_subplot(1, 1, 1)
    ax3 = fig5c.add_subplot(1, 1, 1)
    count = 0
    fs = 14
    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(
                d, basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
                no_flag=True, meta_data=False)
            color = monthly_color[d.month-1]
            smb = symbols[d.month-1]
            s = 'AC240'

            # normal diff
            # scatter = ax1.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].line_amplitude.data - integrated_data['U5303'].line_amplitude.data / , marker=smb, color=color, s=12)
            # ax1.set_ylabel(r'$\Delta T_{B,l}$ [K]', fontsize=fs)

            # frac diff
            scatter = ax1.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].line_amplitude.data - integrated_data['U5303'].line_amplitude.data)/integrated_data['U5303'].line_amplitude.data,  marker=smb, color=color, s=12)
            ax1.set_ylabel(r'$\Delta T_B$ [\%]', fontsize=fs)
            # ax1.set_ylim(-1,0.5)

            ax1.set_title(r'Line amplitude bias $\Delta T_{B,l}$', fontsize=fs+2)
            ax1.set_xlabel(r'Brightness Temperature $T_{B}$ [K]', fontsize=fs)
            theoretical_nonlinearities = np.polyfit(
                [80, 186, 292], [0, -0.20, 0], deg=2)
            fitted_poly_theoretical_nonlinearities = np.poly1d(
                theoretical_nonlinearities)
            ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].continuum_value_line_center.data -
                        integrated_data['U5303'].continuum_value_line_center.data, marker=smb, color=color, s=12)
            ax2.plot(np.arange(70, 300, 1), fitted_poly_theoretical_nonlinearities(
                np.arange(70, 300, 1)), 'k-', linewidth=0.4)
            ax2.axvline(80, color='b', linewidth=0.6, ls='--')
            ax2.axvline(292, color='r', linewidth=0.6, ls='--')
            ax2.set_ylabel(r'$\Delta T_{B,c}$ [K]', fontsize=fs)
            ax2.set_xlabel(r'Brightness Temperature $T_{B}$ [K]', fontsize=fs)
            ax2.set_title(r'Continuum bias $\Delta T_{B,c}$ ', fontsize=fs+2)
            #ax2.set_ylim(-2,0)
            #ax2.set_ylim(-1,1)
            ax3.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].slope_indiv-integrated_data['U5303'].slope_indiv.data)/integrated_data['U5303'].slope_indiv.data,  marker=smb, color=color, s=12)
            ax3.set_ylabel(r'$\Delta m$ [\%]', fontsize=fs)
            #ax3.set_ylim(-0.8,0.2)
            ax3.set_title(r'Slope bias $\Delta m$', fontsize=fs+2)
            ax3.set_xlabel(r'Brightness Temperature $T_B$ [K]', fontsize=fs)
        except:
            print('no data for :', d)
            pass
    legend_elements = [
        Line2D([0], [0], marker='v', color='w', markerfacecolor=monthly_color[0],
               label=month_name[0], markersize=size+2),
        Line2D([0], [0], marker='<', color='w', markerfacecolor=monthly_color[1],
               label=month_name[1], markersize=size+2),
        Line2D([0], [0], marker='^', color='w', markerfacecolor=monthly_color[2],
               label=month_name[2], markersize=size+2),
        Line2D([0], [0], marker='>', color='w', markerfacecolor=monthly_color[3],
               label=month_name[3], markersize=size+2)
    ]
    ax2.text(85, -0.23, '$T_{cold}$', fontsize=fs, color='b')
    ax2.text(262, -0.23, '$T_{hot}$', fontsize=fs, color='r')
    # legend = ax1.legend(*scatter.legend_elements(prop='colors'), month_name, fontsize='xx-small',loc=1, title='Month')
    # ax1.add_artist(legend)
    # ax2.legend(['U5303','AC240'], fontsize='small')
    # ax2.text(85, -0.23, '$T_{cold}$', fontsize=12, color='b')
    # ax2.text(250, -0.23, '$T_{hot}$', fontsize=12, color='r')
    for ax in [ax1,ax2,ax3]:
        ax.grid()
        #ax.set_xlim(80,)
        ax.tick_params(axis='both', which='major', labelsize=14)

    ax1.legend(handles=legend_elements, fontsize=12, loc='lower left')

    #fig2.suptitle('Difference : '+s+' - U5303')
    # fig.suptitle('Mean hot counts')
    fig5a.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig5b.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig5c.tight_layout(rect=[0, 0.01, 1, 0.95])

    plt.show()
    figures2.append(fig5a)
    figures2.append(fig5b)
    figures2.append(fig5c)
    save_single_pdf('/home/eric/Documents/PhD/MOPI/Data/Level3/'+'bias_all_rel' +
                    integration_strategy+'_2021.pdf', figures2)
if plot_bias_TOD_full:
    figures2 = list()
    fig4 = plt.figure(figsize=(9, 6))
    ax1 = fig4.add_subplot(1, 2, 1)
    # ax2 = fig4.add_subplot(1,3,2)
    ax3 = fig4.add_subplot(1, 2, 2)
    count = 0
    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(
                pd.DatetimeIndex([d]), basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
                no_flag=True, meta_data=False)
            color = monthly_color[d.month-1]
            s = 'AC240'
            scatter = ax1.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].line_amplitude.data -
                                                                        integrated_data['U5303'].line_amplitude.data)/integrated_data['U5303'].line_amplitude.data, color=color, s=12)
            # scatter = ax1.scatter(integrated_data[s].time_of_day.data, integrated_data[s].THot.data-integrated_data['U5303'].line_amplitude.data, color=color, s=12)
            ax1.set_ylabel(r'$\Delta T_B$ [%]')
            # ax1.set_ylim(-2,0)
            # ax1.set_ylim(-1,0.5)
            ax1.set_title('Line amplitude difference')
            ax1.set_xlabel(r'Mean Brightness Temperature $T_B$ [K]')
            # ax2.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].continuum_value_line_center.data-integrated_data['U5303'].continuum_value_line_center.data)/integrated_data['U5303'].continuum_value_line_center.data, color=color, s=12)
            # ax2.set_ylabel(r'$\Delta T_B$ [%]')
            # ax2.set_xlabel('Mean $T_B$ [K]')
            # ax2.set_title('$\Delta T_B$ continuum')
            # ax2.set_ylim(-2,0)
            # ax2.set_ylim(-1,1)
            ax3.scatter(integrated_data[s].mean_Tb.data, 100*(integrated_data[s].slope_indiv -
                                                              integrated_data['U5303'].slope_indiv)/integrated_data['U5303'].slope_indiv, color=color, s=12)
            ax3.set_ylabel(r'$\Delta m$ [%]')
            # ax3.set_ylim(-0.8,0.2)
            ax3.set_title(r'Slope difference $\Delta m$')
            ax3.set_xlabel(r'Mean Brightness Temperature $T_B$ [K]')
        except:
            print('no data for :', d)
            pass
    ax1.grid()
    # ax2.grid()
    ax3.grid()
    ax1.legend(handles=legend_elements, fontsize=fs, loc=3)
    fig4.suptitle('Fractional difference with U5303')
    # fig.suptitle('Mean hot counts')
    fig4.tight_layout(rect=[0, 0.01, 1, 0.95])
    plt.show()
    figures2.append(fig4)
    fig3 = plt.figure(figsize=(6, 6))
    # ax1 = fig2.add_subplot(1,1,1)
    ax2 = fig3.add_subplot(1, 1, 1)
    # ax3 = fig2.add_subplot(1,3,3)
    count = 0
    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(
                pd.DatetimeIndex([d]), basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(
                no_flag=True, meta_data=False)
            color = monthly_color[d.month-1]
            s = 'AC240'
            # scatter = ax1.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].line_amplitude.data-integrated_data['U5303'].line_amplitude.data, color=color, s=12)
            # #scatter = ax1.scatter(integrated_data[s].time_of_day.data, integrated_data[s].THot.data-integrated_data['U5303'].line_amplitude.data, color=color, s=12)
            # ax1.set_ylabel(r'$\Delta T_B$ [K]')
            # #ax1.set_ylim(-2,0)
            # #ax1.set_ylim(-1,0.5)
            # ax1.set_title('line amplitude difference')
            # ax1.set_xlabel('Mean $T_B$ [K]')
            theoretical_nonlinearities = np.polyfit(
                [80, 186, 292], [0, -0.20, 0], deg=2)
            fitted_poly_theoretical_nonlinearities = np.poly1d(
                theoretical_nonlinearities)
            ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].continuum_value_line_center.data -
                        integrated_data['U5303'].continuum_value_line_center.data, color=color, s=12)
            ax2.plot(np.arange(70, 303, 1), fitted_poly_theoretical_nonlinearities(
                np.arange(70, 303, 1)), 'k-', linewidth=0.4)
            ax2.axvline(80, color='b', linewidth=0.6, ls='--')
            ax2.axvline(292, color='r', linewidth=0.6, ls='--')
            ax2.set_ylabel(r'$\Delta T_B$ [K]')
            ax2.set_xlabel('Mean $T_B$ [K]')
            ax2.set_title(r'$\Delta T_B$ continuum')
            # ax2.set_ylim(-2,0)
            # ax2.set_ylim(-1,1)
            # ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope_indiv*1e9-integrated_data['U5303'].slope_indiv.data*1e9, color=color, s=12)
            # ax3.set_ylabel('$\Delta m$ [K/GHz]')
            # #ax3.set_ylim(-0.8,0.2)
            # ax3.set_title('Slope difference')
            # ax3.set_xlabel('Mean $T_B$ [K]')
        except:
            print('no data for :', d)
            pass
    legend_elements = [
        Line2D([0], [0], marker='.', color='w', markerfacecolor='magenta',
               label=month_name[0], markersize=size+2),
        Line2D([0], [0], marker='.', color='w', markerfacecolor='blue',
               label=month_name[1], markersize=size+2),
        Line2D([0], [0], marker='.', color='w', markerfacecolor='cyan',
               label=month_name[2], markersize=size+2),
        Line2D([0], [0], marker='.', color='w', markerfacecolor='orange',
               label=month_name[3], markersize=size+2)
    ]
    # legend = ax1.legend(*scatter.legend_elements(prop='colors'), month_name, fontsize='xx-small',loc=1, title='Month')
    # ax1.add_artist(legend)
    # ax2.legend(['U5303','AC240'], fontsize='small')
    ax2.text(85, -0.23, '$T_{cold}$', fontsize=14, color='b')
    ax2.text(270, -0.23, '$T_{hot}$', fontsize=14, color='r')
    # ax1.grid()
    ax2.grid()
    # ax3.grid()
    # ax2.legend(handles=legend_elements,fontsize=14,loc=2)
    # fig2.suptitle('Difference : '+s+' - U5303')
    # fig.suptitle('Mean hot counts')
    fig3.tight_layout(rect=[0, 0.01, 1, 0.95])
    plt.show()
    figures2.append(fig3)
    save_single_pdf(basename_lvl2+'full_bias_all_' +
                    integration_strategy+'_2021.pdf', figures2)
if plot_o3:
    spectro_lvl2 = integration.spectrometers
    level2_data = integration.read_level2(
        spectrometers=spectro_lvl2, extra_base='_all')
    outName = 'bias_o3_feb'
    mopi5_library.plot_O3_all_mopi5(level2_data, outName)

if compare_level2_mopi5:
    from matplotlib import cm
    fs = 16
    spectro_lvl2 = integration.spectrometers
    level2_data = integration.read_level2(
        spectrometers=spectro_lvl2, extra_base='fascod_paper')
    level2_ac240_unbiased = integration.read_level2(
        spectrometers=['AC240'], extra_base='fascodunbiased_all')
    level2_data['AC240_unbiased'] = level2_ac240_unbiased['AC240']

    reference_spectro = 'U5303'

    o3_ref = level2_data[reference_spectro].o3_x.isel(o3_lat=0, o3_lon=0).where(level2_data[reference_spectro].o3_mr.isel(o3_lat=0, o3_lon=0)>0.8)
    o3_z = level2_data[reference_spectro].o3_z

    # Corresponds to altitude range: 20-30km, 30-40km, 40-50km, 50-60km, 60-70km
    p_range = [np.arange(7, 11), np.arange(11, 14), np.arange(
        14, 17), np.arange(17, 20), np.arange(20, 23)]

    p_range_full = np.arange(7, 23)

    bias_AC240 = dict()
    relative_bias = dict()

    mean_o3_alt_range = dict()

    for s in ['AC240', 'AC240_unbiased', 'USRP-A']:
    #for s in ['USRP-A']:
        alt_bias = xr.DataArray()
        mean_o3 = xr.DataArray()
        rel_bias = xr.DataArray()

        o3_s1 = level2_data[s].o3_x.isel(o3_lat=0, o3_lon=0).where(level2_data[s].o3_mr.isel(o3_lat=0, o3_lon=0)>0.8)
        diff = o3_s1 - o3_ref
        #alt_bias =1e6*np.mean(np.abs(diff.isel(o3_p=p_range[0])),1)
        alt_bias = 1e6*np.mean(diff.isel(o3_p=p_range[0]), 1)
        mean_o3 = 1e6*np.mean(o3_ref.isel(o3_p=p_range[0]), 1)
        rel_bias = alt_bias/mean_o3

        alt_bias_full = np.mean(1e6*np.mean(diff.isel(o3_p=p_range_full), 1))
        print('Full bias for ',  s)
        print(alt_bias_full)

        #alt_bias= alt_bias.expand_dims(dim='altitude')
        for i in np.arange(1, len(p_range)):
            #alt_bias = xr.concat([alt_bias, 1e6*np.mean(np.abs(diff.isel(o3_p=p_range[i])),1)], dim='altitude_range')
            alt_bias = xr.concat(
                [alt_bias, 1e6*np.mean(diff.isel(o3_p=p_range[i]), 1)], dim='altitude_range')
            mean_o3 = xr.concat(
                [mean_o3, 1e6*np.mean(o3_ref.isel(o3_p=p_range[i]), 1)], dim='altitude_range')
            re = np.mean(diff.isel(o3_p=p_range[i]), 1)/ np.mean(o3_ref.isel(o3_p=p_range[i]), 1)
            rel_bias = xr.concat([rel_bias, re], dim='altitude_range')

        alt_bias.coords['altitude_range'] = [25, 35, 45, 55, 65]
        alt_bias.coords['time'] = np.arange(0, 15) 
        mean_o3.coords['altitude_range'] = [25, 35, 45, 55, 65]
        mean_o3.coords['time'] = np.arange(0, 15)         
        rel_bias.coords['time'] = np.arange(0, 15)       
        rel_bias.coords['altitude_range'] = [25, 35, 45, 55, 65]  
        bias_AC240[s] = alt_bias
        mean_o3_alt_range[s] = mean_o3
        relative_bias[s] = 100*rel_bias

        plt.fill_betweenx(
            [25, 35, 45, 55, 65], 
            alt_bias.mean(dim='time')-alt_bias.std(dim='time'),
            alt_bias.mean(dim='time')+alt_bias.std(dim='time'), 
            alpha=0.5
        )

        print('mean bias of '+s+' compared to '+reference_spectro +' :', alt_bias.mean(dim='time').data)
        print('corresponding to :',100*alt_bias.mean(dim='time').data / mean_o3.mean(dim='time').data, '%')
        print('std dev of bias of '+s+' compared to '+reference_spectro +' :', alt_bias.std(dim='time').data)

    np.savetxt(outfolder + 'o3_bias_USRP.txt',[bias_AC240['USRP-A'].mean(dim='time').data, 100*bias_AC240['USRP-A'].mean(dim='time').data/mean_o3_alt_range['USRP-A'].mean(dim='time').data , bias_AC240['USRP-A'].std(dim='time').data],fmt='%.2f')
    np.savetxt(outfolder + 'o3_bias_AC240.txt',[bias_AC240['AC240'].mean(dim='time').data,100*bias_AC240['AC240'].mean(dim='time').data/mean_o3_alt_range['AC240'].mean(dim='time').data,bias_AC240['AC240'].std(dim='time').data],fmt='%.2f')
    np.savetxt(outfolder + 'o3_bias_AC240corr.txt',[bias_AC240['AC240_unbiased'].mean(dim='time').data,100*bias_AC240['AC240_unbiased'].mean(dim='time').data/mean_o3_alt_range['AC240_unbiased'].mean(dim='time').data, bias_AC240['AC240_unbiased'].std(dim='time').data],fmt='%.2f')
    toplim = 12
    colormap = 'coolwarm'#cmap_crameri#'bwr'
    colormap = cm.get_cmap('coolwarm')
    colormap.set_bad(color='black')
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=True, figsize=(8, 6))
    pl = relative_bias['AC240'].plot(
        ax=axs[0],
        cmap=colormap,
        center=0,
        vmin=-toplim,
        vmax=toplim,
        # linewidth=0,
        # rasterized=True,
        add_colorbar=False
    )
        #cbar_kwargs={"label": r"$\Delta$O$_3$ [ppmv]"}
    
    pl2 = relative_bias['AC240_unbiased'].plot(
        ax=axs[1],
        cmap=colormap,
        center=0,
        vmin=-toplim,
        vmax=toplim,
        # linewidth=0,
        # rasterized=True,
        add_colorbar=False
    )
    
    #cbar_kwargs={"label": r"$\Delta$O$_3$ [ppmv]"}
    
    # pl3 = bias_AC240['USRP-A'].plot(
    #     ax=axs[2],
    #     cmap=colormap,
    #     center=0,
    #     vmin=-toplim,
    #     vmax=toplim,
    #     # linewidth=0,
    #     # rasterized=True,
    #     cbar_kwargs={"label": r"$\Delta$ O3 [ppm]"}
    # )
    #cb = plt.colorbar(pl, orientation="vertical", pad=0.5)
    cbaxes = fig.add_axes([0.87, 0.25, 0.02, 0.5]) 
    cb = plt.colorbar(pl, cax=cbaxes, orientation="vertical", pad=0.01)
    cb.set_label(label=r"$\Delta$O$_3$ [\%]", size=fs, weight='bold')
    cb.ax.tick_params(labelsize=fs)

    for i in [0,1]:
        axs[i].set_xlabel(r'Mean Brightness Temperature [K]',fontsize=fs)
        axs[i].set_ylabel('Altitude [km]',fontsize=fs)
        axs[i].tick_params(axis='both', which='major', labelsize=fs)
        axs[i].set_xticklabels(['',r'$77$',r'$88$',r'$98$',r'$107$',r'$117$',r'$135$',r'$159$',r'$209$'])

    axs[0].set_title('AC240')
    axs[1].set_title('AC240 corrected')
    #axs[1].set_title('USRP-A')
    pl.set_edgecolor('face')

    plt.tight_layout(rect=[0, 0.03, 0.87, 1])
    fig.savefig('/home/eric/Documents/PhD/MOPI/Data/Level3/' + 'o3_diff_all_rel_feb_tblabel.pdf')

if plot_o3_sel:
    spectro_lvl2 = integration.spectrometers
    level2_data = integration.read_level2(
        spectrometers=spectro_lvl2, extra_base='fascod_paper')
    level2_ac240_unbiased = integration.read_level2(
        spectrometers=['AC240'], extra_base='fascodunbiased_all')
    level2_data['AC240_unbiased'] = level2_ac240_unbiased['AC240']
    outName = '/home/eric/Documents/PhD/MOPI/Data/Level3/' + \
        'bias_o3_feb_all_fascod_fix_noise_rel.pdf'
    # mopi5_library.plot_O3_sel_mopi5(
    #     level2_data, 
    #     spectro=['U5303', 'AC240', 'AC240_unbiased'], 
    #     outName=outName)

    integration.plot_o3_retrieval_mopi5(
        level2_data,
        spectrometers=['U5303', 'AC240', 'AC240_unbiased'],
        idx=[0,11], # np.arange(0,15),#[0,11], 
        save_plot=True,
        identifier=identifier_plot,
        lowerBound=lowerBound,
        outname=outName
    )
# %%
if plot_sel_paper:
    cycles = np.arange(0,15)
    spectro_lvl2 = integration.spectrometers
    level2_data = integration.read_level2(
        spectrometers=spectro_lvl2, extra_base='fascod_paper')
    outname = '/home/eric/Documents/PhD/MOPI/Data/Level3/' +'/'+'o3_comp_3on1_'+integration.datestr + '_plot_all'
    #outname = '/home/esauvageat/Documents/MOPI5/Level3'+'/'+'o3_comp_3on1_'+integration.datestr + '_plot_all'
    mopi5_library.plot_O3_3on1_paper(
        level2_data,
        outname,
        spectrometer=spectro_lvl2,
        cycles= [0]
    )

if plot_avks_paper:
    spectro_lvl2 = integration.spectrometers
    level2_data = integration.read_level2(
        spectrometers=spectro_lvl2, extra_base='fascod_paper')
    outname ='/home/eric/Documents/PhD/MOPI/Data/Level3' +'/'+'o3_comp_avks_'+integration.datestr + '_plot_all'
    mopi5_library.plot_O3_3on1_avks_paper(
        level2_data,
        outname,
        spectrometer=spectro_lvl2,
        cycles= [0]
    )
    