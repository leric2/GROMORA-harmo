#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Integration and DataRetrieval classes implementation for GROMOS instrument


"""
from abc import ABC
import os, datetime
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from xarray.core import dataarray

from gromora_integration import Integration
from gromora_retrievals import DataRetrieval
import GROMORA_library as GROMORA_library

def return_bad_channels_gromos(date):
    '''
    to get the bad channels as a function of the date for GROMOS
    
    Parameters
    ----------
    date : datetime object
    
    '''
    
    bad_channels = np.arange(16383,16384)
    return bad_channels

class IntegrationGROMOSFB(Integration):
    '''
    Implementing the Integration class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, 
        date, 
        basename_lvl1, 
        integration_strategy = None,
        integration_time = 1
        ):

        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "GROMOS"

        self.bandwidth = [1e9]
        spectrometers = ["FB"]
        
        self.lo = 1.45875e11

        level1_folder = os.path.join(basename_lvl1, instrument_name)

        #integration_strategy = 'simple'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder)
    
    def return_bad_channels(self, date, spectro):
        return return_bad_channels_gromos(date)

    # def compare_Tb_chunks(self, dim='time', idx=[0], save = False, Tb_corr = False):
    #     figures = list()

    #     for i in idx:
    #         figures.append(GROSOM_library.plot_Tb_chunks(self, self.integrated_dataset, i)) 

    #     if T_corr:
    #         for i in idx:
    #             figures.append(GROSOM_library.plot_Tb_corr_chunks(self, self.integrated_dataset, i)) 
        
    #     if save:
    #         raise NotImplementedError()

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1'):
        '''
        Correction function for the troposphere. 
        
        Invidual correction for each spectrometers specified !
        '''
        return GROMORA_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

class GROMOS_FB_LvL2(DataRetrieval):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, integration_strategy, integration_time, extra_base=''):
        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "GROMOS"

        self.bandwidth = [1e9]
        spectrometers = ["FB"]

        
        self.lo = 1.45875e11
        self.reference_elevation_angle = 90
        
        level1_folder = basename_lvl1 # os.path.join(basename_lvl1, instrument_name)
        level2_folder = basename_lvl2# os.path.join(basename_lvl2, instrument_name)
        
        self.institution = 'University of Bern;UBERN'
        self.affiliation = 'ubern001'
        self.source = 'MWR.O3_UBERN'
        self.longitude = 7.44
        self.latitude = 46.95
        self.altitude = 560
        
        # Can be used for plotting names (GROMOS_AC240_...)
        self.basename_plot_level2 = instrument_name+'_'+spectrometers[0]+'_'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder, level2_folder, extra_base)
    
    #####################################################################################################################
    # Definition of the retrieval_param dictionary for GROMOS FB
    #####################################################################################################################
    def define_retrieval_param(self, retrieval_param):
        """This function fills the dictionary containing all parameters for GROMOS retrievals.

        Args:
            retrieval_param (dict): a dict containing all required parameters, names, folder path etc for running a GROMORA retrievals.
        """
        retrieval_param["obs_freq"] = self.observation_frequency

        ########################################################
        # Sensor related parameter:
        retrieval_param['sensor'] = 'FB'
        retrieval_param['SB_bias'] = 0
        
        retrieval_param['sideband_response'] = 'theory'
        retrieval_param['use_all_channels'] = False
        retrieval_param['window_corrected_spectrum'] = True
        retrieval_param["f_shift"] = 0

        # Frequency grid for the simulation:
        retrieval_param["number_of_freq_points"] = 4*1201
        retrieval_param["irregularity_f_grid"] = 5000

        ########################################################
        # Pressure grids
        # Computed from altitudes values (all values in meters)

        # Altitude of the surface and radiometer
        retrieval_param["surface_altitude"] = 1000
        retrieval_param["observation_altitude"] =  1000   
        
        # Pressure grid for the simulation (FM)
        retrieval_param["z_top_sim_grid"] = 112e3
        retrieval_param["z_bottom_sim_grid"] = 600 #600
        retrieval_param["z_resolution_sim_grid"] = 2e3

        # Pressure grid for the retrievals 
        retrieval_param["retrieval_grid_type"] = 'altitude'
        retrieval_param["z_top_ret_grid"] = 95e3
        retrieval_param["z_bottom_ret_grid"] = 1e3
        retrieval_param["z_resolution_ret_grid"] = 2e3

        # Pressure value for the continuum retrievals
        retrieval_param["retrieval_h2o_grid_type"] = 'pressure'
        # retrieval_param["z_top_ret_grid_h2o"] = 20e3
        # retrieval_param["z_bottom_ret_grid_h2o"] = 1e3
        # retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

        retrieval_param["h2o_pressure"] = [500e2]

        # Offset for the pointing angle -> instrument_class (usually 0)
        retrieval_param['pointing_angle_corr'] = self.correct_pointing(retrieval_param)


        ########################################################
        # Species definition and spectroscopy
        retrieval_param['spectroscopy_type'] = 'XML'
        # 'H2O-PWR98' #H2O-ContMPM93 'H2O-PWR98, H2O' #"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
        retrieval_param['water_vapor_model'] = 'H2O-PWR98'
        retrieval_param['o2_model'] = 'O2-PWR93'  # 'O2-MPM93'
        # 'N2-SelfContMPM93'
        retrieval_param['n2_model'] = 'N2-SelfContStandardType'

        retrieval_param['selected_species'] = ['O3', retrieval_param['water_vapor_model'],
                                               retrieval_param['o2_model'], retrieval_param['n2_model']]

        #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
        # retrieval_param["azimuth_angle"]=32
        line_file = retrieval_param['ARTS_DATA_PATH'] +"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
        #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"
        #line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all_species.par'
        #line_file = '/home/es19m597/Documents/GROMORA/InputsRetrievals/Hitran_f_modified.xml'
        retrieval_param['line_file'] = line_file
        ########################################################
        # Atmosphere for the simulation
        retrieval_param['atm'] = 'ecmwf_cira86'  # fascod  ecmwf_cira86
        # max_diff, simple_stack_corr, simple, max_diff_surf
        retrieval_param['ptz_merge_method'] = 'max_diff'
        retrieval_param['ptz_merge_max_Tdiff'] = 5
        retrieval_param['ecmwf_store_location'] = '/storage/tub/atmosphere/ecmwf/locations/'+self.location #+str(retrieval_param['date'].year)
        #retrieval_param['ecmwf_store_location'] ='/home/eric/Documents/PhD/ECMWF'
        retrieval_param['cira86_path'] = os.path.join(retrieval_param['ARTS_DATA_PATH'], 'planets/Earth/CIRA86/monthly')
        # time interval [h] around which we collect ECMWF profile
        retrieval_param['extra_time_ecmwf'] = 3.5

        ########################################################
        # A priori and covariances
        retrieval_param['o3_apriori'] = 'waccm_monthly'#'waccm_monthly'
        # 'waccm_yearly_scaled'low_alt_ratio
        retrieval_param['h2o_apriori'] = 'ecmwf'  # 'ecmwf' # 'fascod_extended'
        retrieval_param['o3_apriori_covariance'] = 'sinefit_optimized' # 'low_alt_ratio_optimized' #low_alt_ratio_optimized
        retrieval_param['waccm_file'] = retrieval_param["GROMORA_FOLDER"] + '/files/waccm_o3_climatology.nc' #'/storage/tub/atmosphere/WACCM/waccm_o3_climatology.nc'
        retrieval_param["apriori_O3_cov"] = 1e-6  # 1e-6
        retrieval_param["apriori_H2O_stdDev"] = 1  # 6e-4 12e-4 0.5 16e-4

        ########################################################
        # Polyfit retrieval
        retrieval_param['poly_order'] = 2
        retrieval_param['covmat_polyfit_0'] = 0.1
        retrieval_param['covmat_polyfit_1'] = 0.5
        retrieval_param['covmat_polyfit_2'] = 0.5

        ########################################################
        # Type of noise covariance to use
        retrieval_param['noise_covariance']  = 'noise_level'
        # factor to increase the noise variance
        retrieval_param['increased_var_factor'] = 0.3
        
        ########################################################
        # Variable related to channel response of FB
        retrieval_param['channel_width_measured'] = self.channel_width_measured_good()

        ########################################################
        #  Baseline retrievals
        retrieval_param['sinefit_periods'] = self.baseline_period(retrieval_param) # np.array([400e6]) #np.array([319e6]) 160e6 110e6 119e6 113.5, 66e6, 45e6
        retrieval_param['sinefit_covmat'] = len(retrieval_param['sinefit_periods']) * [np.array([0.1, 0.1])] 
        retrieval_param["binned_ch"] = False

        # OEM parameters, see https://atmtools.github.io/arts-docs-2.4/docserver/methods/OEM.html 
        retrieval_param['oem_method'] = 'gn'# 'lm'
        retrieval_param['max_iter'] = 10
        retrieval_param['stop_dx'] = 20.1
        retrieval_param['lm_ga_setting']=[10, 2.0, 2.0, 100, 1, 99]

        return retrieval_param

    def return_bad_channels(self, date, spectro):
        return return_bad_channels_gromos(date)

    
    def baseline_period(self, retrieval_param):
        '''
        Depending on the dates, function to apply the appropriate baseline periods for the GROMOS retrievals

        ''' 
        #if (retrieval_param['date'] >= datetime.date(2009,1,1)) & (retrieval_param['date'] < datetime.date(2015,2,23)):
        #     baseline_periods = np.array([178e6, 240e6, 360e6])
        # elif  (retrieval_param['date'] >= datetime.date(2015,2,23)) & (retrieval_param['date'] < datetime.date(2015,8,31)):
        #     baseline_periods = np.array([140e6, 240e6, 400e6])
        # elif  (retrieval_param['date'] >= datetime.date(2015,8,31)) & (retrieval_param['date'] < datetime.date(2017,1,1)):
        #     baseline_periods = np.array([160e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2017,1,1)) & (retrieval_param['date'] < datetime.date(2018,1,1)):
        #     baseline_periods = np.array([178e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2018,1,1)) & (retrieval_param['date'] < datetime.date(2019,1,1)):
        #     baseline_periods = np.array([135e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2019,1,1)) & (retrieval_param['date'] < datetime.date(2019,3,15)):
        #     baseline_periods = np.array([155e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2019,3,15)) & (retrieval_param['date'] < datetime.date(2022,2,16)):
        #     baseline_periods = np.array([135e6, 178e6, 240e6])
        # elif retrieval_param['date'] > datetime.date(2022,2,16):
        #     baseline_periods = np.array([135e6, 178e6, 240e6])
        #else:
        baseline_periods = np.array([])
            #raise ValueError('No Sinefit implement for FB yet')
            
        return baseline_periods
    
    def correct_pointing(self, retrieval_param):
        if (retrieval_param['date'] >= datetime.date(2019,2,12)) & (retrieval_param['date'] < datetime.date(2019,3,13)):
            return -5
        else:
            return 0
    ############################################################################################################################
    # Grids definition
    ############################################################################################################################
    # Frequency grid:
    def make_f_grid(self, retrieval_param): 
        """Function to create the frequency grid for the retrievals     

        Args:
            retrieval_param (dict): dictionary with all retrieval parameters 

        Returns:
            f_grid (numpy array): the frequency grid to do the retrievals.
        """

        # Creation of a cubic decay in the resolution around the observation frequency
        n_f = retrieval_param["number_of_freq_points"]  # Number of points
        bw = 1.3*retrieval_param["bandwidth"]  # Bandwidth
        x = np.linspace(-1, 1, n_f)
        f_grid = x ** 3 + x / retrieval_param["irregularity_f_grid"]
        f_grid = f_grid * bw / (max(f_grid) - min(f_grid)) + \
            retrieval_param['obs_freq']
        #f_grid = np.linspace(-1e9, 1e9, n_f)+retrieval_param['obs_freq']
        #f_grid = np.linspace(retrieval_param["f_min"]-10, retrieval_param["f_max"]+10, n_f)
        #f_grid = np.concatenate((f_grid,np.arange(148.975e9,149.975e9,1e6)))
        #f_grid = np.arange(141e9,150e9,1e6)

        # An option to plot the frequency grid:
        if retrieval_param["show_f_grid"]:
            fig = plt.figure()
            plt.semilogy(f_grid[1:]/1e9, np.diff(f_grid)/1e3, '.')
            plt.xlim((retrieval_param['obs_freq']-200e6) /
                     1e9, (retrieval_param['obs_freq']+200e6)/1e9)
            # plt.ylim(0,300)
            plt.ylabel(r'$\Delta f$ [kHz]')
            plt.suptitle('Frequency grid spacing')
            plt.show()

        return f_grid

    def plot_level2(self, ac, ds, retrieval_param, title="",figure_list = list()):
        '''
        Plotting function directly taken from Jonas ;)

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
        fshift_ret = None
        if retrieval_param['retrieval_quantities'] == 'o3_h2o':
            ozone_ret, h2o_ret = ac.retrieval_quantities
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift':
            ozone_ret, h2o_ret, fshift_ret = ac.retrieval_quantities
            #print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
            print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_polyfit':
            ozone_ret, h2o_ret, polyfit_ret = ac.retrieval_quantities
            print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit':
            ozone_ret, h2o_ret, polyfit_ret, fshift_ret = ac.retrieval_quantities
            print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
            print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit_sinefit':
            ozone_ret, h2o_ret, polyfit_ret, fshift_ret, sinefit_ret = ac.retrieval_quantities
            print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
            print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
            print('Sinefit : ', sinefit_ret.x)
        else:
            ozone_ret,  = ac.retrieval_quantities

        good_channels = ds.good_channels[retrieval_param['integration_cycle']].data == 1
       # f_backend = ds.frequencies[retrieval_param['integration_cycle']].values[good_channels]
        f_backend = ozone_ret.ws.f_backend.value
       # y = ds.Tb[retrieval_param['integration_cycle']].values[good_channels]
        y = ozone_ret.ws.y.value
        yf = ac.yf[0]
        r = y - yf
        bl = ac.ws.y_baseline.value 

        # Modelling error:
        e_mod = np.matmul(ozone_ret.ws.dxdy.value, r)[0:len(ozone_ret.p_grid)]

        fig, axs = plt.subplots(2, sharex=True, figsize=(12,10))
        axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, y, label="observed")
        axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, yf, label="fitted")
        #axs[0].set_ylim(-5, 50)
        axs[0].legend()
        if fshift_ret is not None:
            axs[0].text(
                0.02,
                0.8,
                'fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3),
                transform=axs[0].transAxes,
                verticalalignment="bottom",
                horizontalalignment="left",
        )

        axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
        axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, bl, label="baseline")
        axs[1].set_ylim(np.median(r)-2.5, np.median(r+2.5))
        axs[1].legend()
        axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(retrieval_param['obs_freq'] / 1e9))

        for ax in axs:
            ax.set_ylabel("$T_B$ [K]")
            ax.set_xlim([min((f_backend - retrieval_param['obs_freq']) / 1e6), max((f_backend - retrieval_param['obs_freq']) / 1e6)])
        fig.suptitle(title)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        figure_list.append(fig)

        fig, axs = plt.subplots(1, 2, sharey=True, figsize=(12,10))
        axs[0].plot(
            ozone_ret.x * 1e6, ozone_ret.p_grid / 1e2, label="retrieved", marker="x"
        )
        axs[0].plot(ozone_ret.xa * 1e6, ozone_ret.p_grid / 1e2, label="apriori")
        axs[0].set_xlim(-0.1,8.5)
        axs[0].invert_yaxis()
        axs[0].set_yscale('log')
        axs[0].set_xlabel("Ozone VMR [ppm]")
        axs[0].set_ylabel("Pressure [hPa]")
        axs[0].legend()
        axs[1].plot(ozone_ret.mr, ozone_ret.p_grid / 1e2)
        axs[1].set_xlabel("Measurement response")
        #axs[1].invert_yaxis()
        #axs[1].set_yscale('log')

        axs[0].grid(True)
        axs[1].grid(True)
        figure_list.append(fig)

        fig, axs = plt.subplots(1, 2, sharey=True, figsize=(12,10))    

        # Modelling error:
        error_mod = np.matmul(ozone_ret.ws.dxdy.value, r)[0:len(ozone_ret.p_grid)]

        axs[0].plot(ozone_ret.es * 1e6, ozone_ret.p_grid / 1e2, label="smoothing error")
        axs[0].plot(ozone_ret.eo * 1e6, ozone_ret.p_grid / 1e2, label="obs error")
        axs[0].plot(error_mod * 1e6, ozone_ret.p_grid / 1e2, label="modelling error")

        axs[0].set_xlabel("$e$ [ppm]")
        axs[0].set_ylabel("Pressure [hPa]")
        axs[0].legend()
        axs[0].set_xlim(-0.5,1)

        # axs[1].plot(100*(ozone_ret.x - og_ozone)/og_ozone, ozone_ret.z_grid / 1e3, label="retrieval-og")
        # axs[1].plot(100*(ozone_ret.xa - og_ozone)/og_ozone, z_og / 1e3, label="apriori-og")
        # axs[0].set_xlabel("Rel diff [%]")
        # axs[0].set_ylabel("Altitude [km]")

        for avk in ozone_ret.avkm:
            if 0.8 <= np.sum(avk) <= 1.2:
                axs[1].plot(avk, ozone_ret.p_grid / 1e2)

        axs[1].set_xlabel("AVKM")
        axs[0].invert_yaxis()
        axs[0].set_yscale('log')

        #axs[1].invert_yaxis()
        #axs[1].set_yscale('log')

        axs[1].grid(True)
        axs[0].grid(True)

        fig.suptitle(title + " Ozone")
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        figure_list.append(fig)

        #print('Atmospheric opacity :')
        #print(ac.ws.y_aux.value, ', mean: ', np.mean(ac.ws.y_aux.value))

        if retrieval_param['plot_opacities']:
            #opacities = ds.tropospheric_opacity[retrieval_param['integration_cycle']].values[good_channels]
            #plt.plot(f_backend,opacities,label='matlab' )
            fig, ax = plt.subplots(1, 1, sharey=True)
            ax.plot(f_backend,ac.ws.y_aux.value[0], label='ARTS' )
            ax.legend()

        if (retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift') or (retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit'):
            h2o_apriori = ac.ws.vmr_field.value[1,:,0,0]
            h2o_profile = h2o_ret.x[0]*ac.ws.vmr_field.value[1,:,0,0]
            z_h2o = ac.ws.z_field.value[:,0,0]

            fig, axs = plt.subplots(2, 2, sharey=True, figsize=(12,10))
            axs[0][0].semilogx(
                h2o_profile, z_h2o / 1e3, label="retrieved", marker="x"
            )
            axs[0][0].semilogx(h2o_apriori, z_h2o / 1e3, label="apriori")
            axs[0][0].set_xlabel("Water vapor [rel]")
            axs[0][0].set_ylabel("Altitude [km]")
            axs[0][0].legend()

            axs[0][1].plot(h2o_ret.mr, h2o_ret.z_grid / 1e3, marker="x")
            axs[0][1].set_xlabel("Measurement response")

            axs[1][0].plot(h2o_ret.es, h2o_ret.z_grid / 1e3, label="smoothing error", marker="x")
            axs[1][0].plot(h2o_ret.eo, h2o_ret.z_grid / 1e3, label="obs error", marker="x")
            axs[1][0].set_xlabel("$e$ [rel]")
            axs[1][0].set_ylabel("Altitude [km]")
            axs[1][0].legend()

            for avk in h2o_ret.avkm:
                #if 0.8 <= np.sum(avk) <= 1.2:
                axs[1][1].plot(avk, h2o_ret.z_grid / 1e3)
            axs[1][1].set_xlabel("AVKM")

            axs[1][1].set_ylim((0,20))
            axs[0][0].grid(True)
            axs[0][1].grid(True)
            axs[1][1].grid(True)
            axs[1][0].grid(True)

            # axs[0][0].set_ylim(-0.5, 30)
            # axs[0][1].set_ylim(-0.5, 30)
            # axs[1][1].set_ylim(-0.5, 30)
            # axs[1][0].set_ylim(-0.5, 30)

        #axs[0][0].grid(True)
        #axs[0][1].grid(True)
        #axs[1][1].grid(True)
        #axs[1][0].grid(True)

            fig.suptitle(" Water vapor retrieval (v{})".format(1))
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            figure_list.append(fig)

            temp = ac.ws.t_field_raw.value.to_xarray()
            alt = ac.ws.z_field_raw.value.to_xarray()
            fig, axs = plt.subplots(1, 2, figsize=(12,10))
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
    
            figure_list.append(fig)   

        return figure_list

    # @property
    # def day2flag_level2(self):
    #     '''
    #     A selection of days to flags for the level2 GROMOS data. 
    #     These days have been identified in the GROMORA time series detailed analysis that can be found in the GROMORA retrievals UG.

    #     '''
    #     date2flag_gromos =  [
    #         datetime.date(2015,8,26), datetime.date(2015,8,27), datetime.date(2015,8,28),
    #         pd.date_range('2012-07-24', '2012-08-07'),
    #         pd.date_range('2019-01-14', '2019-02-12'),

    #     ]
    #     return date2flag_gromos

    def cost_threshold(self, year):
        '''
        Cost threshold over which we flag the level 2
        '''
        return 5

    def channel_width_measured(self):
        '''
        Measured by Axel
        '''
        return 1e6*np.array([100, 114.667784,118.600345,109.055542,146.597255,164.5163,18.906618,18.189048,21.780333,10.595859,5.630094,6.079148,2.830508,2.815202,3.103099,1.794136,1.39955,0.826034,0.770904,0.242118,0.230111,0.225406,0.231984,0.243707,0.23766,0.231579,0.235805,0.233243,0.646687,0.642863,1.560209,1.261859,3.162788,2.869584,3.370188,6.127845,7.779886,11.394174,22.580517,23.252718,31.629541,178.777324,170.339139,180.648311,139.313079,159.515409,180.338123,151.604634])
    
    def channel_width_measured_good(self):
        '''
        TO DO, clean it for bad channels
        '''
        good_channel_width = self.channel_width_measured()
    
        return good_channel_width

    def channel_width_theory(self):
        return 1e6*[100,100,100,100,100,100,30,20,20,10,5,5,2,2,2,1,1,0.5,0.5,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,1,1,2,2,2,5,5,10,20,20,30,100,100,100,100,100,100,100]

    @property
    def usb_grid(self):
        return np.arange(148.975e9,150.175e9,100e6)

    @property
    def basecolor(self):
       return '#d7191c' 

    @property
    def polyfit_threshold(self):
        return 0.1

    @property
    def standard_air_pressure(self):
        return 955
        
    @property
    def standard_air_temperature(self):
        return 10

    @property
    def cycle_duration(self):
        return 7/3600
    
    @property
    def global_attributes_ndacc(self):
        pi_mail='axel.murk@unibe.ch'
        do_name = 'Sauvageat;Eric'
        do_mail = 'eric.sauvageat@unibe.ch'
        rou= 'Please contact Axel Murk at axel.murk@unibe.ch'
        ackn ='The middle atmospheric ozone radiometer GROMOS is operated by the Institute of Applied Physics, University of Bern, Switzerland.'
        description='Atmospheric ozone profiles from continuous measurements of ground-based 142 GHz-microwave radiometer GROMOS at Bern, Switzerland'
        return dict(
            PI_NAME=self.name_PI,
            PI_AFFILIATION=self.institution,
            PI_ADDRESS = self.contact,
            PI_EMAIL = pi_mail,
            DO_NAME= do_name,
            DO_AFFILIATION=self.institution,
            DO_ADDRESS = self.contact,
            DO_EMAIL =  do_mail,
            DS_NAME= do_name,
            DS_AFFILIATION=self.institution,
            DS_ADDRESS = self.contact,
            DS_EMAIL = do_mail,
            DATA_DESCRIPTION = description,
            DATA_RULES_OF_USE =rou,
            DATA_ACKNOWLEDGEMENT=ackn,
            FILE_PROJECT_ID = 'NDACC-GROMOS'
        )