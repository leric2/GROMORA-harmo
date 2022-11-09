#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 

@author: eric

Collection of functions for dealing with time in the frame of GROMORA retrievals

"""
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
from pytz import timezone, utc

# local time zone for GROMOS and SOMORA:
local_timezone = timezone('Europe/Zurich')

# time zone of the measurement times from GROMOS and SOMORA:
gromora_tz = timezone('UTC')


def solar_zenith_angle(ha, doy, lat):
    # sun declination
    declination = -23.44*np.cos(np.deg2rad((doy+10)*360/365))
    cos_declination = np.cos(np.deg2rad(declination))
    sin_declination = np.sin(np.deg2rad(declination))

    # Hour angle 
    cos_solar_hour_angle = np.cos(np.deg2rad(ha))

    # Sunrise/Sunset:
    cos_hour_angle_night = -np.tan(np.deg2rad(lat))*np.tan(np.deg2rad(declination))

    # if cos_solar_hour_angle < cos_hour_angle_night:
    #     night = True
    # else:
    #     night = False
    cos_sza = sin_declination*np.sin(np.deg2rad(lat)) + np.cos(np.deg2rad(lat))*cos_declination*cos_solar_hour_angle
    sza = np.rad2deg(np.arccos(cos_sza))

    return sza

def hour_angle_sunset(doy, lat):
    declination = -23.44*np.cos(np.deg2rad((doy+10)*360/365))
    cos_declination = np.cos(np.deg2rad(declination))
    sin_declination = np.sin(np.deg2rad(declination))
    cos_lat = np.cos(np.deg2rad(lat))
    # Sunrise/Sunset:
    cos_hour_angle_sunset = -np.tan(np.deg2rad(lat))*np.tan(np.deg2rad(declination))    

    # NOAA formula:
    cos_hour_angle_sunrise_sunset = np.cos(np.deg2rad(90.833))/(cos_lat*cos_declination) - np.tan(np.deg2rad(lat))*np.tan(np.deg2rad(declination))

    hour_angle_sunrise = np.rad2deg(np.arccos(cos_hour_angle_sunset))

    hour_angle_sunrise_NOAA = np.rad2deg(np.arccos(cos_hour_angle_sunrise_sunset))

    return hour_angle_sunrise, hour_angle_sunrise_NOAA

def pysolar_sza(date, lat, lon):
    # Using pysolar package:
    time = utc.localize(datetime64_2_datetime(date))
    sza_pysolar = solar.get_altitude(lat, lon, time)
    return sza_pysolar

def datetime64_2_datetime(dt64):
    unix_epoch = np.datetime64(0, 's')
    one_second = np.timedelta64(1, 's')
    seconds_since_epoch = (dt64 - unix_epoch) / one_second
    date = dt.datetime.utcfromtimestamp(seconds_since_epoch)
    return date

def equation_of_time(doy):
    B=(360/365)*(doy-81)
    eot = 9.87*np.sin(np.deg2rad(2*B)) - 7.53*np.cos(np.deg2rad(B)) -1.5*np.sin(np.deg2rad(B))
    return eot

def time_correction_factor(lon, lstm, eot):
    return 4*(lon - lstm) + eot

def local_solar_time(local_time, tc):
    return local_time + tc/60

def hour_angle(lst):
    hours_from_midnight = (lst - lst.replace(hour=0, minute=0,second=0, microsecond=0)).total_seconds()/3600
    return 15*(hours_from_midnight-12)

def lst_sunset_from_hour_angle(sunset_ha, midnight_lst):
    sunrise = midnight_lst + dt.timedelta(hours=12) + dt.timedelta(hours=-np.abs(sunset_ha)/15)
    sunset = midnight_lst + dt.timedelta(hours=12) + dt.timedelta(hours=np.abs(sunset_ha)/15)
    return sunrise, sunset

def get_sunset_lst_from_lst(lst, lat):
    sunset_ha, sunset_NOAA = hour_angle_sunset(lst,lat)

    sunrise_lst, sunset_lst = lst_sunset_from_hour_angle(
        sunset_ha, 
        midnight_lst = lst.replace(hour=0, minute=0,second=0, microsecond=0)
    )
    return sunrise_lst, sunset_lst

def get_LST_from_GROMORA(date, lat, lon, print_option=False, check_format=True):
    #dt = utc.localize(datetime64_2_datetime(date))
    if check_format:
        if np.issubdtype(date.dtype, np.datetime64):
            date = datetime64_2_datetime(date).replace(tzinfo=gromora_tz)

    if print_option:
        print('UTC time: ',date) 
    local_time =  date.astimezone(local_timezone)
    if print_option:
        print('Local time: ',local_time)

    doy = pd.to_datetime(date).dayofyear

    eot = equation_of_time(doy)
    if print_option:
        print('Equation of time : ', str(eot))

    lstm = 15*local_time.utcoffset().seconds/3600
    tc = time_correction_factor(lon, lstm, eot)

    lst = local_time + dt.timedelta(minutes=tc)
    
    lst = lst.replace(tzinfo=None)

    ha = hour_angle(lst)

    ha_sunset, ha_NOAA= hour_angle_sunset(doy, lat)
    if print_option:
        print('Hour angle: ', str(ha))
    #  print('Hour angle sunset: ', str(ha_sunset))
    # print('Hour angle NOAA: ', str(ha_NOAA))
    
    if np.abs(ha) > np.abs(ha_NOAA):
        night = True
    else:
        night = False

    sza = solar_zenith_angle(ha, doy, lat)
    return lst, ha, sza, night, tc

def mjd2k_date(dates):
    """
    % ABSTRACT  | Function to compute the MJD2K time. Taken from GEOMS report 1.0:
    %           | https://avdc.gsfc.nasa.gov/PDF/GEOMS/geoms-1.0.pdf
    %           |
    %           |
    %           |
    % ARGUMENTS | INPUTS: a pandas datetime timestamp
    %           |
    %           | OUTPUTS: - mjd2k_date: a double number corresponding to the
    %           | Modified Julian Datetime 2000. It corresponds to the number
    %           | of day from 01.01.2000. 
    """
    mjd2k_date = list()
    for date in dates:
        month = date.month

        if month>2:
            y = date.year
            m = month-3
            d = date.day
        else:
            y = date.year-1
            m = month+9
            d = date.day

        j = int(365.25*(y+4712)) + int(30.6*m+0.5) + 58.5 + d
        if j < 2299159.5:
            jd = j
        else:
            gn = 38 - int(3*int(49 + y/100)/4)
            jd = j+gn

        df = (date.hour*3600.0 + date.minute*60.0 + date.second)/86400
        mjd2k_date.append(jd + df - 2451544.5)
    
    return mjd2k_date


if __name__ == "__main__":
    #date = spectro_dataset.time.data[0]
    date = np.datetime64('2021-10-20 00:06:00.123456')
    
    dt = datetime64_2_datetime(date).replace(tzinfo=gromora_tz)
    doy= pd.to_datetime(date).dayofyear

    lat = np.array(46.95)
    lon = np.array(7.44)

    sza_pysolar_elevation_angle = pysolar_sza(date, lat, lon)

    lst, ha, sza, night, tc = get_LST_from_GROMORA(dt, lat, lon)
    lst_simone = (dt + dt.timedelta(hours=lon*24/360)).replace(tzinfo=None)

    LT_from_lst = (lst - dt.timedelta(minutes=tc))

    solar_noon = (lst.replace(hour=12, minute=0,second=0, microsecond=0) - dt.timedelta(minutes=tc))

    print('Local solar time: ',lst)
    print('Local solar time (Simone): ',lst_simone)
    sza_simone = solar_zenith_angle(hour_angle(lst_simone), doy, lat)

    print('Solar zenith angle (pysolar) = ', 90-sza_pysolar_elevation_angle)
    print('Solar zenith angle (corrected) = ', sza)
    print('Solar zenith angle (Simone) = ', sza_simone)

    print('LST diff:', (lst-lst_simone.replace(tzinfo=None)).total_seconds()/60, 'min')
    if night:
        print('Night !')
    else:
        print('Day !')

    sunset_ha, sunset_NOAA = hour_angle_sunset(doy, lat)

    sunrise_lst, sunset_lst = lst_sunset_from_hour_angle(
        sunset_ha, 
        midnight_lst = lst.replace(hour=0, minute=0,second=0, microsecond=0)
    )
    print('Approx. sunrise (LST): ',sunrise_lst, ' and sunset (LST): ', sunset_lst)
    print('Approx. sunrise (LT): ',(sunrise_lst - dt.timedelta(minutes=tc)), ' and sunset (LT): ', (sunset_lst- dt.timedelta(minutes=tc)))

    print('Daylight hours: ',(sunset_lst-sunrise_lst).total_seconds()/3600, 'hours')
    if np.abs(ha) > sunset_ha:
        print('Night !')
    else:
        print('Day !')