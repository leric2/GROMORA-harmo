#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

This file contains the old definition of the DataRetrieval Abstract class.

"""
from abc import ABC
import os
import datetime

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from gromora_utils import save_single_pdf

from dotenv import load_dotenv

# For ARTS, we need to specify some paths
#load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

import GROMORA_library, retrieval_module


class DataRetrieval(ABC):
    def __init__(
        self,
        instrument_name=None,
        observation_frequency=None,
        spectrometers=None,
        integration_strategy='classic',
        integration_time=None,
        date=None,
        level1_folder=None,
        level2_folder=None,
        extra_base=None
        ):
        
        self.instrument_name = instrument_name
        self.date = date
        #self.datestr = date.strftime('%Y_%m_%d')
        self.observation_frequency = observation_frequency
        self.spectrometers = spectrometers
        self.integration_strategy = integration_strategy
        self.int_time = integration_time
        self.level1_folder = level1_folder
        self.level2_folder = level2_folder

        # must be false for Integration
        self.multiday = False
        try:
            if len(self.date) > 1:
                self.multiday = True
        except TypeError:
            self.multiday = False
            
    
        self.calibrated_data = dict()
        self.meteo_data = dict()
        self.calibration_flags = dict()
        self.filename_level1a = dict()
        self.filename_level1b = dict()
        #extra_base =''
        for s in self.spectrometers:
            self.filename_level1a[s] = list()
            self.filename_level1b[s] = list()
            if self.multiday:
                for d in self.date:
                    self.datestr = d.strftime('%Y_%m_%d')
                    self.filename_level1a[s].append(
                        os.path.join(
                    self.level1_folder,
                    self.instrument_name + "_level1a_" +
                    s + "_" + self.datestr + extra_base
                    ))
                    if self.integration_strategy == 'classic':
                        if self.int_time == 1:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_" +
                            s + "_" + self.datestr + extra_base
                            )
                        else:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                            s + "_" + self.datestr + extra_base
                            )
                    else:
                        self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                            s + "_" + self.datestr + extra_base
                            )
            else:
                self.datestr = self.date.strftime('%Y_%m_%d')
                if self.integration_strategy == 'classic':
                    if self.int_time == 1:
                        self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_" +
                        s + "_" + self.datestr + extra_base
                        ))
                    else:
                        self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                        s + "_" + self.datestr + extra_base
                        ))
                else:
                    self.datestr = self.date.strftime('%Y_%m_%d')
                    self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                        s + "_" + self.datestr + extra_base
                        ))

    def get_hot_load_temperature(self, time):
        ''' Get hot load temperature for a specific time.
        '''
        raise NotImplementedError("abstract base class")

    def read_level1b(self, no_flag=False, meta_data=True, extra_base=None):
        ''' 
        Reading level1b dataset and completing the information on the instrument
        
        '''
        self.integrated_data = dict()
        self.flags = dict()
        self.integrated_meteo = dict()
        #self.filename_level1a = dict()
        #self.filename_level1b = dict()
        self.filename_level2 = dict()
        for i, s in enumerate(self.spectrometers):
            self.filename_level1a[s] = os.path.join(
            self.level1_folder,
            self.instrument_name + "_level1a_" +
            s + "_" + self.datestr
            )
            if self.int_time == 1:
                # self.filename_level1b[s] = os.path.join(
                # self.level1_folder,
                # self.instrument_name + "_level1b_" +
                # s + "_" + self.datestr
                # )

                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                s + "_" + self.datestr + extra_base
                )
            else:
                # self.filename_level1b[s] = os.path.join(
                # self.level1_folder,
                # self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                # s + "_" + self.datestr
                # )

                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" + str(self.int_time) +"h_" +
                s + "_" + self.datestr + extra_base
                )
        
            #print('reading : ', self.filename_level1b[s])
            #self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = GROSOM_library.read_level1(self._filename_lvl1b)
            try:
                self.integrated_data[s], self.flags[s], self.integrated_meteo[s], global_attrs_level1b = GROSOM_library.read_level1(self.filename_level1b[s][0], no_flag=no_flag)
                print('Read : ', self.filename_level1b[s][0])
            except FileNotFoundError:
                print('No file for this day, skipping ', self.filename_level1b[s][0])


        if meta_data:       
            # Meta data
            self.institution = global_attrs_level1b['institution']
            self.instrument_name_from_level1b = global_attrs_level1b['instrument']
            self.location = global_attrs_level1b['location']

            self.number_of_spectrometer = global_attrs_level1b['number_of_spectrometer']
            self.calibration_version = global_attrs_level1b['calibration_version']

            #self.raw_data_filename = global_attrs_level1b['raw_data_filename']
            self.raw_data_software_version = global_attrs_level1b['raw_data_software_version']

            self.filename_level1a = global_attrs_level1b['filename_level1a']
            self.raw_file_warning = global_attrs_level1b['raw_file_warning']
            self.labview_logfile_warning = global_attrs_level1b['labview_logfile_warning']

        # some information from the ds
        #self.number_of_channels = len(self.level1b_ds.channel_idx.values)

        #self.frequencies = self.level1b_ds.frequencies.values
        #self.IF = self.level1b_ds.intermediate_freq.values
        if self.multiday:
            print('No meta data updates for multi-days reading (only first day is saved) !')
            for i in range(1,len(self.date)):
                for s in self.spectrometers:
                    try:
                        calibrated_data, calibration_flags, meteo_data, global_attrs_level1b = GROSOM_library.read_level1(self.filename_level1b[s][i], no_flag=no_flag)
                    except FileNotFoundError:
                        print('No file for this day, skipping ', self.filename_level1b[s][i])
                    else:
                        print('Read : ', self.filename_level1b[s][i])
                        self.calibrated_data[s] = xr.concat([self.calibrated_data[s],calibrated_data], dim='time')
                        self.calibration_flags[s] = xr.concat([self.calibration_flags[s],calibration_flags], dim='time')
                        self.meteo_data[s] = xr.concat([self.meteo_data[s], meteo_data], dim='time')
                        
        #self.time = self.level1b_ds.time.values
        #self.number_of_time_records = len(self.time)
        return self.integrated_data, self.flags, self.integrated_meteo
        #return self.level1b_ds, self.flags, self.meteo_ds

    def read_level2(self, spectrometers=None, extra_base=''):
        ''' 
        Reading level2 dataset and completing the information on the instrument
        
        '''
        if spectrometers is None:
            spectrometers = self.spectrometers

        if self.multiday:
            self.integrated_data = dict()
            self.filename_level2 = dict()
            self.level2_data=dict()
            print('No meta data updates for multi-days reading (only first day is saved) !')
            counter = 0
            for d in self.date:
                for s in spectrometers:
                    if self.int_time == 1:
                        self.filename_level2[s] = os.path.join(
                        self.level2_folder,
                        self.instrument_name + "_level2_" +
                        s + "_" + d.strftime('%Y_%m_%d') + extra_base
                        )
                    else:
                        self.filename_level2[s] = os.path.join(
                        self.level2_folder,
                        self.instrument_name + "_level2_" + str(self.int_time)+'h_' +
                        s + "_" + d.strftime('%Y_%m_%d') + extra_base
                        )
                    print(self.filename_level2[s]) 
                    try:
                        level2_data = xr.open_dataset(
                            self.filename_level2[s] + ".nc",
                        )                    
                    except FileNotFoundError:
                        print('No file for this day, skipping ', self.filename_level2[s])
                       # level2_data[s] = xr.Dataset()
                    else:
                        if counter == 0:
                            self.level2_data[s] = level2_data
                        else:
                            self.level2_data[s] = xr.concat([self.level2_data[s],level2_data], dim='time')
                        print('Read : ', self.filename_level2[s])
                        counter = counter + 1
        else:
            self.integrated_data = dict()
            self.filename_level2 = dict()
            self.level2_data=dict()
            for s in spectrometers:
                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                s + "_" + self.datestr + extra_base
                )

                try:
                    self.level2_data[s] = xr.open_dataset(
                        self.filename_level2[s] + ".nc",
                        )
                    print('Read : ', self.filename_level2[s])
                except FileNotFoundError:
                    print('No file for this day, skipping ', self.filename_level2[s])
        return self.level2_data

    def add_bias_characterization_with_non_linearities(self, spectrometers, dim, param_slope, around_center_value):
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

        # For each spectrometer individually
        for s in spectrometers:
            theoretical_nonlinearities =  np.polyfit([80,186,292],[0,-0.20,0],deg=2)
            fitted_poly_theoretical_nonlinearities = np.poly1d(theoretical_nonlinearities) 
        
            non_lin=fitted_poly_theoretical_nonlinearities(self.integrated_data[s].Tb.data)
            
            #meanTb_lc = self.integrated_data[s].Tb.where(abs(self.integrated_data[s].bin_freq_interp-self.observation_frequency)<around_center_value,drop=True).mean(dim='bin_freq_interp')
            # bias_lc = meanTb_lc-meanTb_lc_basis
            # bias_lc_fract = bias_lc/meanTb_lc_basis

            # #self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc)
            # self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc_fract)

            cleanTb = (self.integrated_data[s].Tb+fitted_poly_theoretical_nonlinearities(self.integrated_data[s].Tb)).where(self.integrated_data[s].good_channels==1)
            # Tb_diff = cleanTb-cleanTb_basis
            # Tb_diff_fract = Tb_diff/cleanTb_basis
            right_wing_center = param_slope[s][0]
            right_wing_interval = param_slope[s][1]
            left_wing_center = param_slope[s][2]
            left_wing_interval = param_slope[s][3]

            mean_Tb_right_wing = cleanTb.where(abs(self.integrated_data[s].frequencies[0]-right_wing_center)<right_wing_interval,drop=True).mean(dim='channel_idx')
            mean_Tb_left_wing = cleanTb.where(abs(self.integrated_data[s].frequencies[0]-left_wing_center)<left_wing_interval,drop=True).mean(dim='channel_idx')

            #mean_diff_right_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-right_wing_center)<right_wing_interval,drop=True).mean(dim='bin_freq_interp')
            #mean_diff_left_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-left_wing_center)<left_wing_interval,drop=True).mean(dim='bin_freq_interp')

            slope = (mean_Tb_right_wing - mean_Tb_left_wing) / (right_wing_center - left_wing_center)

            self.integrated_data[s] = self.integrated_data[s].assign(slope_indiv_nonlin = slope)
            
            self.integrated_data[s].slope_indiv.attrs['left_f']=left_wing_center
            self.integrated_data[s].slope_indiv.attrs['right_f']=right_wing_center
            
            Tb0 = mean_Tb_left_wing-slope*left_wing_center

            continuum_value_line_center = self.observation_frequency * slope + Tb0

            self.integrated_data[s] = self.integrated_data[s].assign(continuum_value_line_center_nonlin = continuum_value_line_center)

            line_center_amplitude = cleanTb.where(abs(self.integrated_data[s].frequencies-self.observation_frequency)<around_center_value,drop=True).mean(dim='channel_idx')
            diff_lc_continuum = line_center_amplitude - continuum_value_line_center
            self.integrated_data[s] = self.integrated_data[s].assign(line_amplitude_nonlin = diff_lc_continuum)

        return self.integrated_data

    def plot_level1b_TB(self, title='', save=False, outfolder='', save_name='int_spectra', idx=None):
        figures = list()
        
        if idx is None:
            figures.append(GROSOM_library.plot_Tb_all(self, self.integrated_data, title=title)) 
        else:
            figures.append(GROSOM_library.plot_Tb_selected(self, self.integrated_data, title=title, idx=idx)) 

        if save:
            save_single_pdf(outfolder+self.instrument_name+'/'+save_name+self.datestr+'.pdf', figures)
            print('saved in '+outfolder+self.instrument_name+'/'+save_name+self.datestr+'.pdf')
            #save_pngs(self.level1_folder+'time_series_'+self.datestr+'_', figures)


    # def plot_level1b_TB(self, calibration_cycle):
    #     plt.plot(self.integrated_data[self.spectrometers].frequencies,level1b_dataset.Tb_trop_corr[calibration_cycle])
    #     plt.ylim((0,200))
    #     pass
    
    def plot_meteo_ds_level1b_dataset(self):
        GROSOM_library.plot_meteo_level1b(self.meteo_ds)
        pass 
    
    def forward_model(self, retrieval_param):
        ''' 
        Performing a FM according to the retrieval_param dictionnary as input
        parameters
        
        '''
        return retrieval_module.forward_model(retrieval_param)
    
    def retrieve_cycle(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None, ac=None, sensor = None):
        ''' 
        Performing single retrieval for a given calibration cycle (defined in retrieval_param) 
        '''
        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False

        retrieval_param['ref_elevation_angle'] =  90

        return retrieval_module.retrieve_cycle(self, spectro_dataset, retrieval_param, ac_FM=ac, sensor=sensor)

    def retrieve_daily(self, spectro_dataset, retrieval_param):
        ''' 
        Trying parallel retrievals for the cycles defined in retrieval_param
        '''
        retrieval_param["binned_ch"] = False

        retrieval_param['ref_elevation_angle'] =  90

        return retrieval_module.retrieve_daily(self, spectro_dataset, retrieval_param)

    def retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None, ac=None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False
 
        retrieval_param['ref_elevation_angle'] = 90

        return retrieval_module.retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, ac_FM=ac)

    def test_retrieval(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        retrieval_param["binned_ch"] = False
 
        retrieval_param['ref_elevation_angle'] = 90

        return retrieval_module.test_retrieval(spectro_dataset, retrieval_param)

    def smooth_and_apply_correction(self, level1b_dataset, meteo_ds):   
        '''
        doing quick tropospheric correction as it was not saved in lvl1b

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        meteo_ds : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        return GROSOM_library.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    
    def create_binning(self, freq, tb, retrieval_param):   
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        meteo_ds : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        n_f = retrieval_param["number_of_freq_points"]
        bw_extra = 1.2*self.bandwidth

        return GROSOM_library.create_bin_vector(self.observation_frequency, freq, tb, n_f, bw_extra)

    def bin_spectrum(self, freq, tb, bin_vect):   
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        meteo_ds : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''

        return GROSOM_library.bin_spectrum(freq, tb, bin_vect)

    def smooth_corr_spectra(self, level1b_dataset, retrieval_param):
        '''
        
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        return GROSOM_library.smooth_corr_spectra(level1b_dataset, retrieval_param)
    
    def find_bad_channels(self, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        return GROSOM_library.find_bad_channels(self.level1b_ds, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
    
    def find_bad_channels_stdTb(self, spectrometers, stdTb_threshold, apply_on='cal', dimension=['time','channel_idx']):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if apply_on=='cal':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.date, spectro)
                self.calibrated_data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.calibrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.calibrated_data
        elif apply_on=='int':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.date, spectro)
                self.integrated_data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.integrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.integrated_data
        else:
            raise ValueError()

    def find_bad_channels_stdTb_old(self, spectro, stdTb_threshold):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        bad_channels = self.return_bad_channels(date, spectro)
        self.data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.data[spectro], bad_channels, stdTb_threshold)

        return self.data[spectro]
    
    def plot_level2_from_tropospheric_corrected_spectra(self, ac, spectro_dataset, retrieval_param, title, figure_list):
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        ac : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.
        title : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        return GROSOM_library.plot_level2_from_tropospheric_corrected_spectra(spectro_dataset, ac, retrieval_param, title, figure_list)
    
    def plot_level2(self, ac, spectro_dataset, retrieval_param, title, figure_list):
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        ac : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.
        title : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        return GROSOM_library.plot_level2(spectro_dataset, ac, retrieval_param, title, figure_list)

    def plot_level2_test_retrieval(self, ac, retrieval_param, title):
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        ac : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.
        title : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        return GROSOM_library.plot_level2_test_retrieval(ac, retrieval_param, title)
    
    def return_bad_channels(self, date, spectro):
        '''
        it is overwritten by each instrument classes
        '''
        return []

def run_retrieval(instrument,retrieval_param):
    a=2
    '''
    In this function we call the retrieval process step-by-step
    
    instrument is actually a concrete implementation of the the abstract 
    DataRetrieval class --> will be initiated here in the future
    
    In the end, we will store this function elsewhere along with the classes definition.
    '''
    return a

if __name__=="__main__":
    
    instrument_name = "GROMOS"
    date = datetime.date(2019,2,1)
    int_time = 1

    basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    
    #basename_lvl1 = "/scratch/GROSOM/Level1/"
    #basename_lvl2 = "/scratch/GROSOM/Level2/"
    
    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
    #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"

    if instrument_name=="GROMOS":
        from gromos_classes import IntegrationGROMOS, GROMOS_LvL2
        calibration = IntegrationGROMOS(date, basename_lvl1, basename_lvl2, int_time)
        instrument = GROMOS_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    elif instrument_name=="SOMORA":
        from somora_classes import SOMORA_LvL2
        instrument = SOMORA_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    elif instrument_name=="mopi5":
        from mopi5_classes import MOPI5_LvL2
        #basename_lvl1 = "/scratch/MOPI5/Level1/"
        #basename_lvl2 = "/scratch/MOPI5/Level2/"
        basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
        basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        instrument = MOPI5_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    
    # Dictionnary containing all EXTERNAL retrieval parameters 
    retrieval_param = dict()

    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param["retrieval_type"] = 1

    retrieval_param["obs_freq"] = instrument.observation_frequency
    
    retrieval_param["plot_meteo_ds"] = True
    retrieval_param["number_of_freq_points"] = 601
    retrieval_param["irregularity_f_grid"] = 10
    retrieval_param["show_f_grid"] = True

    retrieval_param["z_top_sim_grid"] = 97e3
    retrieval_param["z_bottom_sim_grid"] = 800
    retrieval_param["z_resolution_sim_grid"] = 1e3

    retrieval_param["z_top_ret_grid"] = 95e3
    retrieval_param["z_bottom_ret_grid"] = 800
    retrieval_param["z_resolution_ret_grid"] = 3e3

    retrieval_param["z_top_ret_grid_h2o"] = 50e3
    retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

    retrieval_param['increased_var_factor'] = 1
    retrieval_param['unit_var_y']  = 3


    retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
    retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    retrieval_param["apriori_O3_cov"] = 1.5e-6

    retrieval_param['water_vapor_model'] = "H2O-PWR98"
    #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
    #retrieval_param["azimuth_angle"]=32
    
    #retrieval_param['obs_freq'] = 1.4217504e11
    retrieval_param['line_file'] = line_file
    retrieval_param['atm'] ='fascod'
    retrieval_param['ecmwf_store_location'] ='/scratch/ECMWF'
    retrieval_param['extra_time_ecmwf'] = 6

    retrieval_param['cira86_path'] = os.path.join(
        ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')
    # Check the structure of the file and maybe use it ?
    #print(netCDF4.Dataset(filename+".nc").groups.keys())
    
    # level1b_dataset, meteo_ds, global_attr_lvl1b_ds = run_retrieval(instrument,retrievalParam)
    calibrated_data, cal_flags, cal_meteo = calibration.read_level1a()
    data, flags, meteo = instrument.read_level1b()
    assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'

    if instrument_name == 'mopi5':
        instrument = instrument.correction_function_mopi5('U5303', 290)
        instrument.plot_comparison_mopi5_spectrometers(calibration_cycle=[0,1,2,3])
         
    #for i,s in enumerate(instrument.spectrometer):
    spectro = 'AC240'
    spectro_dataset = instrument.data[spectro]
    #retrieval_param = {**global_attrs_level1b, **retrieval_param}
    #else :
    #    raise ValueError('incoherent instrument definition')
    
    #level1b_dataset = instrument.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    #instrument.return_bad_channels(date,'U5303')
    #spectro_dataset = instrument.find_bad_channels_stdTb(spectro, 20)

    # instrument.data = instrument.data.find_bad_channels(
    #     spectrometer=spectro,
    #     Tb_min=0,
    #     Tb_max=260,
    #     boxcar_size=128,
    #     boxcar_thresh=7
    # )
    
    # instrument.plot_meteo_ds_level1b_dataset()

    #level1b_dataset = instrument.smooth_corr_spectra(level1b_dataset, retrieval_param)
    #f_sim, y_sim = instrument.forward_model(retrieval_param)
    #plt.plot(f_sim, y_sim[0], level1b_dataset.frequencies.values, level1b_dataset.Tb[1].values)
    
    retrieval_param["integration_cycle"] = 0

    #bin_vector = instrument.create_binning(
    #    instrument.frequencies,
    #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
    #    retrieval_param
    #)

    #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)

    
    if retrieval_param["retrieval_type"] == 1:
        retrieval_param["surface_altitude"] = 10e3
        retrieval_param["observation_altitude"] =  15e3
        ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
        figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(ac, spectro_dataset, retrieval_param, title = 'retrieval_trop_corr')
        level2 = ac.get_level2_xarray()
        #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
        save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
    elif retrieval_param["retrieval_type"] == 2:
        retrieval_param["surface_altitude"] = 1200
        retrieval_param["observation_altitude"] =  1200
        ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
        figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='retrieval_o3_h20')
        save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_'+retrieval_param['water_vapor_model']+'.pdf', figure_list)
    elif retrieval_param["retrieval_type"] == 3:
        retrieval_param["surface_altitude"] = 1500
        retrieval_param["observation_altitude"] =  1500
        ac, retrieval_param = instrument.test_retrieval(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
        figure_list = instrument.plot_level2_test_retrieval(ac, retrieval_param, title ='test_retrieval_o3_h2o')
        save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'Perrin_with_h2o.pdf', figure_list)
    else:
        pass