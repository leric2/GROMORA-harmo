#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:28:46 2020

@author: eric

Collection of functions for dealing with GROSOM data

Including : 
    * level1a, level1b and level2 data
    * a-priori data
"""

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

#from retrievals import arts

def read_level1(filenameLevel1):
    """Example function with types documented in the docstring.
    Description HERE

    Args:
        param1 (int): The first parameter.
        param2 (str): The second parameter.

    Returns:
        bool: The return value. True for success, False otherwise.

    """
    
    DS = xr.open_dataset(
        filenameLevel1 + ".nc",
        group="spectrometer1",
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        #use_cftime=True,
        )

    globalAttributes=xr.open_dataset(filenameLevel1 + ".nc").attrs
    
    METEO=xr.open_dataset(
        filenameLevel1+".nc",
        group="meteo",
        decode_times=True,
        decode_coords=True,
        )

    flags=xr.open_dataset(
        filenameLevel1+".nc",
        group="flags",
        decode_times=True,
        decode_coords=True,
        )
        
    return DS, flags, METEO, globalAttributes

def correct_troposphere(calibration, spectrometers, dim, method='Ingold_v1'):
    '''
    generic tropospheric correction for GROSOM
    '''
    # little trick to work with other dimensions than time
    
    for s in spectrometers:
        if dim != 'time':
            calibration.integrated_dataset[s] = calibration.integrated_dataset[s].swap_dims({dim:'time'})
        
        tb_corr_da = xr.DataArray(None,
            coords = [calibration.integrated_dataset[s].coords['time'], calibration.integrated_dataset[s].coords['channel_idx']])
        mean_opacity = xr.DataArray(None,
            coords = [calibration.integrated_dataset[s].coords['time']])
        mean_transmittance = xr.DataArray(None,
            coords = [calibration.integrated_dataset[s].coords['time']])
        for i, t in enumerate(calibration.integrated_dataset[s].coords['time']):
            try:
                temperature = calibration.integrated_meteo[s].air_temperature.sel(time=t).data.item()
            except:
                print('No ground temperature found, taking default')
                temperature = 293
            if method == 'Ingold_v1':
                delta_T = 10.4
                tropospheric_temperature = temperature - delta_T
            elif method == 'Ingold_v2':
                raise NotImplementedError
            
            bad_channels = calibration.integrated_dataset[s].good_channels.sel(time=t).data
            bad_channels[bad_channels==0]=np.nan
            clean_tb = calibration.integrated_dataset[s].Tb.sel(time=t).data * bad_channels
            
            corr_func = tropospheric_correction(
                f = calibration.integrated_dataset[s].frequencies.sel(time=t).data,
                y = clean_tb,
                t_trop=tropospheric_temperature,
                use_wings='both',
                skip=[100,100],
                num_el=500)
            
            tb_corr, opac, transm = corr_func(
                f = calibration.integrated_dataset[s].frequencies.sel(time=t).data,
                y = calibration.integrated_dataset[s].Tb.sel(time=t).data
                )
            tb_corr_da.loc[dict(time=t)] = tb_corr
            mean_transmittance.loc[dict(time=t)] = transm
            mean_opacity.loc[dict(time=t)] = opac
        calibration.integrated_dataset[s] = calibration.integrated_dataset[s].assign(Tb_corr = tb_corr_da)
            
        calibration.integrated_dataset[s] = calibration.integrated_dataset[s].assign(tropospheric_transmittance = mean_transmittance)
            
        calibration.integrated_dataset[s] = calibration.integrated_dataset[s].assign(tropospheric_opacity = mean_opacity)
        
        if dim != 'time':
            calibration.integrated_dataset[s] = calibration.integrated_dataset[s].swap_dims({'time':dim})
    return calibration.integrated_dataset


def tropospheric_correction(f, y, t_trop, use_wings='both', skip=[100,100], num_el=500):
    """
    Generate a correction function, given frequencies, brightness temperatures
    and a weighted mean torpospheric temperature.
    Returns a function f: (f, y) -> y_corr

    modified from the OG one from Jonas
    """
    #f = spectro_dataset.frequencies[cycle].data
    #y = spectro_dataset.Tb[cycle].data

    # Skip skip[0] channels at beginning and skip[1] in end
    # use 500 channels for reference
    wings = (slice(skip[0], skip[0]+num_el), slice(-(skip[1]+num_el), -skip[1]))

    #yw = np.stack([y[s].mean() for s in wings])
    #fw = np.stack([f[s].mean() for s in wings])
    yw = np.stack([np.nanmean(y[s]) for s in wings])
    fw = np.stack([np.nanmean(f[s]) for s in wings])

    # Fit a line with slope to the spectrum
    # y = a*y0 + b
    a = (yw[1] - yw[0]) / (fw[1] - fw[0])
    b = yw[1] - a * fw[1]

    #print(a)
    #print(b)

    # Correct by frequency dependant opacity
    def correct(f, y):
        y_eff = a * f + b

        transmittance = (t_trop - y_eff) / (t_trop - 2.7)
        opacity = -np.log(transmittance)
        y_corr = (y - t_trop * (1 - transmittance)) / transmittance
        return y_corr, np.nanmean(opacity), np.nanmean(transmittance)

    #return y_corr, np.nanmean(opacity), np.nanmean(transmittance)
    return correct

def find_bad_channels(level1b_dataset, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh):
    '''
    daily processing
    '''
    good_channels = np.zeros((len(level1b_dataset.time),len(level1b_dataset.channel_idx)))
        
    # identify additional spurious channels on this day    
    for i in range(len(level1b_dataset.time)):
        values = level1b_dataset.Tb[i].values
        smoothed_val = np.convolve(values, np.ones((boxcar_size,)) / boxcar_size, mode="same")
        
        if sum(np.isnan(values))>0:
            print('Some NANs are present in cycle {}'.format(i))
            values[np.isnan(values)] = -9999
            smoothed_val[np.isnan(smoothed_val)] = -9999
        
        boxcar_filter_val = np.abs(values - smoothed_val)
        good_channels_i = np.all([(values < Tb_max),(values > Tb_min),(boxcar_filter_val < boxcar_thresh)], axis=0)
               
        good_channels[i,:] = good_channels_i
        
    # Some known spurious channels (all spectra)
    
    good_channels[:,bad_channels] = 0
    
    
    level1b_dataset = level1b_dataset.assign(
        good_channels = xr.DataArray(good_channels, dims = ['time', 'channel_idx']))
    
    return level1b_dataset

def find_bad_channels_stdTb(level1b_dataset, bad_channels, stdTb_threshold = 20, dimension=['time','channel_idx']):
    '''
    daily processing
    '''
    #good_channels = np.zeros((len(level1b_dataset.time),len(level1b_dataset.channel_idx)))
    good_channels = np.zeros((level1b_dataset.dims[dimension[0]], level1b_dataset.dims[dimension[1]]))
        
    # identify additional spurious channels on this day    
    for i in range(level1b_dataset.dims[dimension[0]]):
        good_channels_i = 1*(level1b_dataset.stdTb[i].values < stdTb_threshold)
        
        good_channels[i,:] = good_channels_i
        
    # Some known spurious channels (all spectra)
    if bad_channels.size != 0:
        good_channels[:,bad_channels] = 0
    
    
    level1b_dataset = level1b_dataset.assign(
        good_channels = xr.DataArray(good_channels, dims = dimension))
    
    return level1b_dataset

def smooth_corr_spectra(level1b_dataset, retrieval_param):
    '''
    

    Parameters
    ----------
    level1b_dataset : TYPE
        DESCRIPTION.
    retrieval_param : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    boxcar_size = retrieval_param['boxcar_size']
    
    # good_indices = np.ones(len(level1b_dataset.Tb[1]))
    
    Tb_trop_smoothed = np.ones((len(level1b_dataset.time),len(level1b_dataset.Tb[1])))*np.nan
    for i in range(len(level1b_dataset.time)):
        r = level1b_dataset.Tb_trop_corr[i].values
        Tb_trop_smoothed[i,:] = np.convolve(r, np.ones((boxcar_size,)) / boxcar_size, mode="same")
    
    level1b_dataset = level1b_dataset.assign(
        Tb_trop_corr_smoothed = xr.DataArray(Tb_trop_smoothed, dims = ['time', 'channel_idx']))
    
    return level1b_dataset

def smooth_and_apply_correction(level1b_dataset, meteo_ds):
    '''
    Function coming from pywdp, in "troposphere.py" and adapted to GROSOM

    Parameters
    ----------
    
    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    '''
    #Ttrop = get_tropospheric_temperature(ds, mean_temperature_func)
    tau = level1b_dataset.trospheric_opacity.values
    T_trop = meteo_ds.air_temperature.values
    
    Tb_tropospheric_corr_smoothed = np.ones((len(level1b_dataset.time),len(level1b_dataset.Tb[1])))*np.nan
    Tb_tropospheric_corr = np.ones((len(level1b_dataset.time),len(level1b_dataset.Tb[1])))*np.nan
    for i in range(len(level1b_dataset.time)):
        Tb_smoothed = np.convolve(level1b_dataset.Tb[i].values, np.ones((1035,)) / 1035, mode="same")
        Tb_tropospheric_corr_smoothed[i,] = (Tb_smoothed - T_trop[i] * (1 - xr.ufuncs.exp(-tau[i]))) / xr.ufuncs.exp(-tau[i])
        Tb_tropospheric_corr[i,:] = (level1b_dataset.Tb[i].values - T_trop[i] * (1 - xr.ufuncs.exp(-tau[i]))) / xr.ufuncs.exp(-tau[i])
    
    tbco = xr.DataArray(Tb_tropospheric_corr, dims = ['time', 'channel_idx'])
    tbco_smoothed = xr.DataArray(Tb_tropospheric_corr_smoothed, dims = ['time', 'channel_idx'])
    
    level1b_dataset = level1b_dataset.assign(Tb_trop_corr = tbco)
    level1b_dataset = level1b_dataset.assign(Tb_trop_corr_smoothed = tbco_smoothed)
    return level1b_dataset

def bin_spectrum(f ,tb ,bin_vect):
    newFreq=np.ones(len(bin_vect)) * np.nan
    newTb = np.ones(len(bin_vect)) * np.nan

    ch_down = 0
    for i in range(len(bin_vect)):
        upCh = int(bin_vect[i])
        if ch_down == upCh:
            newFreq[i] = f[ch_down]
            newTb[i] = tb[ch_down]
        else:     
            newFreq[i] = np.nanmean(f[ch_down:upCh])
            newTb[i] = np.nanmean(tb[ch_down:upCh])
        ch_down = upCh+1
    
    return newFreq, newTb

def create_bin_vector(F0, f, tb, n_f, bw_extra):    
    n_ch = len(f)

    df = (max(f)-min(f))/n_ch

    x = np.linspace(-1, 1, n_f)

    f_grid = (x ** 3 + x/250)
    f_grid = f_grid*bw_extra / (max(f_grid) - min(f_grid)) + F0
    #f_grid = f_grid[np.where(np.logical_and(f_grid>min(f), f_grid<max(f)))]

    diff_f_grid = np.diff(f_grid)

    newFreq = np.ones(len(f_grid)) * np.nan
    newTb = np.ones(len(f_grid)) * np.nan
    binned_ch = np.ones(len(f_grid)) * np.nan
    bin_vect = np.ones(len(f_grid)) * np.nan

    for i in range(len(f_grid)):
        newF = f_grid[i]
        if i ==0:
            channels_2_bin = np.where(f<newF+0.5*diff_f_grid[i])
        elif i==len(f_grid)-1:
            channels_2_bin = np.where(f>=newF-0.5*diff_f_grid[i-1])
        else:
            channels_2_bin = np.where(np.logical_and(f>=newF-0.5*diff_f_grid[i-1], f<newF+0.5*diff_f_grid[i]))
        if len(channels_2_bin[0]) > 0:
            newTb[i] = np.mean(tb[channels_2_bin])
            binned_ch[i] = len(channels_2_bin[0])
            newFreq[i] = np.mean(f[channels_2_bin])
            bin_vect[i] = channels_2_bin[0][-1]
    
    newFreq=newFreq[~np.isnan(binned_ch)]
    newTb=newTb[~np.isnan(binned_ch)]
    bin_vect=bin_vect[~np.isnan(binned_ch)]
    binned_ch=binned_ch[~np.isnan(binned_ch)]
    
    fig, axs = plt.subplots(3)
    axs[0].plot(f/1e9,tb)
    axs[0].plot(newFreq/1e9,newTb)
    axs[0].set_ylim(60,105)
    axs[1].plot(newFreq/1e9,binned_ch)
    axs[1].set_ylim(-1,200)
    axs[2].plot(newFreq/1e9,binned_ch)
    axs[2].set_xlim((F0-50e6)/1e9,(F0+50e6)/1e9)
    axs[2].set_ylim(0,35)
    
    return bin_vect

def plot_Tb_chunks(self, ds_dict, calibration_cycle):
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

def plot_Tb_corr_chunks(self, ds_dict, calibration_cycle):
    fig, axs = plt.subplots(1,1,sharex=True)
    for s in self.spectrometers:
        mask = ds_dict[s].good_channels[calibration_cycle].data
        mask[mask==0]=np.nan
        axs.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,ds_dict[s].Tb_corr[calibration_cycle]*mask, lw=0.5, label=s)
        #axs.set_xlim(110.25, 111.4)
        axs.set_ylim(-5,40)
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
    
def plot_meteo_level1b(METEO):
    """Example function with types documented in the docstring.
    Description HERE

    Args:
        param1 (int): The first parameter.
        param2 (str): The second parameter.

    Returns:
        bool: The return value. True for success, False otherwise.

    """
    fig, axs = plt.subplots(3, sharex=True)
    axs[0].plot(METEO.time.values,METEO.air_temperature.values-273.15) 
    axs[0].set_ylabel('$T_{air}$ [°C]')
    axs[0].set_ylim(-10, 30)
    
    axs[1].plot(METEO.time.values,METEO.relative_humidity.values) 
    axs[1].set_ylabel('rel_hum [-]')
    axs[1].set_ylim(0, 1)
    
    axs[2].plot(METEO.time.values,METEO.precipitation.values) 
    axs[2].set_ylabel('rain [mm]')
    #axs[2].set_ylim(0, 5)
    
    axs[2].set_xlabel('time [h]')
    
    fig.suptitle("Meteo Data")

def plot_level2(ds, ac, retrieval_param, title=""):
    '''
    Plotting function directly taken from Jonas ;)
    OG can be found in retrieval.py in MOPI retrievals
    
    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    ac : TYPE
        DESCRIPTION.
    title : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    figures : TYPE
        DESCRIPTION.

    '''
    
    ozone_ret, h2o_ret = ac.retrieval_quantities
    good_channels = ds.good_channels[retrieval_param['integration_cycle']].data == 1
    f_backend = ds.frequencies.values[good_channels]
    y = ds.Tb[retrieval_param['integration_cycle']].values[good_channels]
    #y = ac.y[0]
    yf = ac.yf[0]
    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")
    
    figures = list()

    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, y, label="observed")
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, yf, label="fitted")
    #axs[0].set_ylim(-5, 50)
    axs[0].legend()
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r_smooth, label="residuals smooth")
    #axs[1].set_ylim(-2, 2)
    axs[1].legend()
    axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(retrieval_param['obs_freq'] / 1e9))

    for ax in axs:
        ax.set_ylabel("$T_B$ [K]")
        ax.set_xlim([min((f_backend - retrieval_param['obs_freq']) / 1e6), max((f_backend - retrieval_param['obs_freq']) / 1e6)])
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].plot(
        ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlim(-0.5,15)
    axs[0][0].set_xlabel("Ozone VMR [ppm]")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")
    
    axs[1][0].plot(ozone_ret.es * 1e6, ozone_ret.z_grid / 1e3, label="smoothing error")
    axs[1][0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3, label="obs error")
    axs[1][0].set_xlabel("$e$ [ppm]")
    axs[1][0].set_ylabel("Altitude [km]")
    axs[1][0].legend()

    for avk in ozone_ret.avkm:
        #if 0.8 <= np.sum(avk) <= 1.2:
        axs[1][1].plot(avk, ozone_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")

    axs[0][0].grid(True)
    axs[0][1].grid(True)
    axs[1][1].grid(True)
    axs[1][0].grid(True)

    fig.suptitle(" Ozone retrieval")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)
    
    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].semilogx(
        h2o_ret.x*1e6, h2o_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].semilogx(h2o_ret.xa*1e6, h2o_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlabel("Water VMR [ppm]")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(h2o_ret.mr, h2o_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")

    axs[1][0].semilogx(h2o_ret.es*1e6, h2o_ret.z_grid / 1e3, label="smoothing error")
    axs[1][0].semilogx(h2o_ret.eo*1e6, h2o_ret.z_grid / 1e3, label="obs error")
    axs[1][0].set_xlabel("$e$ [ppm]")
    axs[1][0].set_ylabel("Altitude [km]")
    axs[1][0].legend()

    for avk in h2o_ret.avkm:
        #if 0.8 <= np.sum(avk) <= 1.2:
        axs[1][1].plot(avk, h2o_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")


    axs[0][0].grid(True)
    axs[0][1].grid(True)
    axs[1][1].grid(True)
    axs[1][0].grid(True)

    axs[0][0].set_ylim(-0.5, 30)
    axs[0][1].set_ylim(-0.5, 30)
    axs[1][1].set_ylim(-0.5, 30)
    axs[1][0].set_ylim(-0.5, 30)
        
    #axs[0][0].grid(True)
    #axs[0][1].grid(True)
    #axs[1][1].grid(True)
    #axs[1][0].grid(True)

    fig.suptitle(" Water vapor retrieval (v{})".format(1))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    temp = ac.ws.t_field_raw.value.to_xarray()
    alt = ac.ws.z_field_raw.value.to_xarray()
    fig, axs = plt.subplots(1, 2)
    axs[0].plot(temp.sel(Latitude=0, Longitude=0).data, temp.Pressure)
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('T [K]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(temp.sel(Latitude=0, Longitude=0).data, alt.sel(Latitude=0, Longitude=0).data/1e3)
    #axs[1].invert_yaxis()
    #axs[1].set_yscale('log')
    axs[1].set_xlabel('T [K]')
    axs[1].set_ylabel('$Z$ [km]')
    fig.suptitle('Raw PTZ profile')
   
    figures.append(fig)   

    
    return figures

def plot_level2_from_tropospheric_corrected(ds, ac, retrieval_param, title=""):
    '''
    Plotting function directly taken from Jonas ;)
    OG can be found in retrieval.py in MOPI retrievals
    
    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    ac : TYPE
        DESCRIPTION.
    title : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    figures : TYPE
        DESCRIPTION.

    '''
    ozone_ret = ac.retrieval_quantities[0]
    good_channels = ds.good_channels[retrieval_param['integration_cycle']].data == 1
    f_backend = ds.frequencies.values[good_channels]
    y = ds.Tb_corr[retrieval_param['integration_cycle']].values[good_channels]
    #y = ac.y[0]
    yf = ac.yf[0]
    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")

    figures = list()

    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, y, label="observed")
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, yf, label="fitted")
    axs[0].set_ylim(-5, 50)
    axs[0].legend()
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r_smooth, label="residuals smooth")
    axs[1].set_ylim(-2, 2)
    axs[1].legend()
    axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(retrieval_param['obs_freq'] / 1e9))

    for ax in axs:
        ax.set_ylabel("$T_B$ [K]")
        ax.set_xlim([min((f_backend - retrieval_param['obs_freq']) / 1e6), max((f_backend - retrieval_param['obs_freq']) / 1e6)])
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].plot(
        ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlim(-0.5,15)
    axs[0][0].set_xlabel("Ozone VMR [ppm]")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")
    
    axs[1][0].plot(ozone_ret.es * 1e6, ozone_ret.z_grid / 1e3, label="smoothing error")
    axs[1][0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3, label="obs error")
    axs[1][0].set_xlabel("$e$ [ppm]")
    axs[1][0].set_ylabel("Altitude [km]")
    axs[1][0].legend()

    for avk in ozone_ret.avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1][1].plot(avk, ozone_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")

    axs[0][0].grid(True)
    axs[0][1].grid(True)
    axs[1][1].grid(True)
    axs[1][0].grid(True)


    fig.suptitle(" Ozone")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)


    temp = ac.ws.t_field_raw.value.to_xarray()
    alt = ac.ws.z_field_raw.value.to_xarray()
    fig, axs = plt.subplots(1, 2)
    axs[0].plot(temp.sel(Latitude=0, Longitude=0).data, temp.Pressure)
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('T [K]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(temp.sel(Latitude=0, Longitude=0).data, alt.sel(Latitude=0, Longitude=0).data/1e3)
    #axs[1].invert_yaxis()
    #axs[1].set_yscale('log')
    axs[1].set_xlabel('T [K]')
    axs[1].set_ylabel('$Z$ [km]')
    fig.suptitle('Raw PTZ profile')
   
    figures.append(fig)   

    return figures

def plot_level2_test_retrieval(ac, retrieval_param, title=""):
    '''
    Plotting function directly taken from Jonas ;)
    OG can be found in retrieval.py in MOPI retrievals
    
    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    ac : TYPE
        DESCRIPTION.
    title : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    figures : TYPE
        DESCRIPTION.

    '''
    ozone_ret, h2o_ret = ac.retrieval_quantities
    #good_channels = ds.good_channels[retrieval_param['integration_cycle']].data == 1
    f_backend = ac.ws.y_f.value
    y = ac.y[0]
    #y = ac.y[0]
    yf = ac.yf[0]
    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")

    figures = list()

    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, y, label="observed")
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, yf, label="fitted")
    axs[0].set_xlim(-0.5,15)
    axs[0].legend()
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r_smooth, label="residuals smooth")
    #axs[1].set_ylim(-2, 2)
    axs[1].legend()
    axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(retrieval_param['obs_freq'] / 1e9))

    for ax in axs:
        ax.set_ylabel("$T_B$ [K]")
        ax.set_xlim([min((f_backend - retrieval_param['obs_freq']) / 1e6), max((f_backend - retrieval_param['obs_freq']) / 1e6)])
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].plot(
        ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlabel("Ozone VMR [ppm]")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")
    
    axs[1][0].plot(ozone_ret.es * 1e6, ozone_ret.z_grid / 1e3, label="smoothing error")
    axs[1][0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3, label="obs error")
    axs[1][0].set_xlabel("$e$ [ppm]")
    axs[1][0].set_ylabel("Altitude [km]")
    axs[1][0].legend()

    for avk in ozone_ret.avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1][1].plot(avk, ozone_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")

    axs[0][0].grid(True)
    axs[0][1].grid(True)
    axs[1][1].grid(True)
    axs[1][0].grid(True)


    fig.suptitle(title + " Ozone")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)


    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].semilogx(
        h2o_ret.x, h2o_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].semilogx(h2o_ret.xa, h2o_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlabel("Water VMR []")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(h2o_ret.mr, h2o_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")

    axs[1][0].plot(h2o_ret.es, h2o_ret.z_grid / 1e3, label="smoothing error")
    axs[1][0].plot(h2o_ret.eo, h2o_ret.z_grid / 1e3, label="obs error")
    axs[1][0].set_xlabel("$e$ []")
    axs[1][0].set_ylabel("Altitude [km]")
    axs[1][0].legend()

    for avk in h2o_ret.avkm:
        #if 0.8 <= np.sum(avk) <= 1.2:
        axs[1][1].plot(avk, h2o_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")

    axs[0][0].grid(True)
    axs[0][1].grid(True)
    axs[1][1].grid(True)
    axs[1][0].grid(True)

    fig.suptitle(title + " Water (v{})".format(1))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    temp = ac.ws.t_field_raw.value.to_xarray()
    alt = ac.ws.z_field_raw.value.to_xarray()
    fig, axs = plt.subplots(1, 2)
    axs[0].plot(temp.sel(Latitude=0, Longitude=0).data, temp.Pressure)
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('T [K]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(temp.sel(Latitude=0, Longitude=0).data, alt.sel(Latitude=0, Longitude=0).data/1e3)
    #axs[1].invert_yaxis()
    #axs[1].set_yscale('log')
    axs[1].set_xlabel('T [K]')
    axs[1].set_ylabel('$Z$ [km]')
    fig.suptitle('Raw apriori ptz profile')
   
    figures.append(fig)   

    
    return figures
    