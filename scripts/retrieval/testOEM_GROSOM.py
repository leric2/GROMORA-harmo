"""
This is very close to the "controlfiles/artscomponents/TestOEM.arts" cfile shipped with arts.
Plus some plots to show the retrieval results.

Inspired from Jonas, this is a script to test the OEM retrievals for GROMOS.

Ozone retrieval only on tropospheric corrected spectra.
"""
import retrieval_module
import apriori_data_GROSOM
import data_GROSOM
from retrievals import covmat
from retrievals import arts
from retrievals.arts.atmosphere import p2z_simple, z2p_simple
from typhon.arts.workspace import arts_agenda
from dotenv import load_dotenv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore', message='numpy.dtype size changed')

# load_dotenv('/home/esauvageat/Documents/ARTS/arts-examples/.env')
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

show_plots = True
save_plots = True
save_netcdf = True

name = 'testOEM_GROMOS_ECMWF_apriori'

# For testing
basename = "/home/eric/Documents/PhD/GROSOM/Level1/"
level2_data_folder = "/home/eric/Documents/PhD/GROSOM/Level2/"

# basename="/home/esauvageat/Documents/GROSOM/Analysis/Level1/GROMOS/"
# level2_data_folder = "/home/esauvageat/Documents/GROSOM/Analysis/Level2/"

line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
# line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"

'''
Define to make the retrieval work as if launched from main_test.py,
therefore using the retrieval_param structure even if not really needed...
'''
@arts_agenda
def inversion_iterate_agenda(ws):
    """Custom inversion iterate agenda to ignore bad partition functions."""
    ws.Ignore(ws.inversion_iteration_counter)

    # Map x to ARTS' variables
    ws.x2artsAtmAndSurf()
    ws.x2artsSensor()

    # To be safe, rerun some checks
    ws.atmfields_checkedCalc(negative_vmr_ok=0)
    ws.atmgeom_checkedCalc()

    # Calculate yf and Jacobian matching x
    ws.yCalc(y=ws.yf)

    # Add baseline term
    ws.VectorAddVector(ws.yf, ws.y, ws.y_baseline)

    # This method takes cares of some "fixes" that are needed to get the Jacobian
    # right for iterative solutions. No need to call this WSM for linear inversions.
    ws.jacobianAdjustAndTransform()

instrument_name = "GROMOS"
filename = basename+"GROMOS_level1b_AC240_2019_04_16"

retrieval_param = dict()
retrieval_param["integration_cycle"] = 6
retrieval_param["plot_meteo_ds"] = True
retrieval_param["number_of_freq_points"] = 601


retrieval_param["azimuth_angle"] = 32
retrieval_param["observation_altitude"] = 15e3
retrieval_param['obs_freq'] = 1.4217504e11
retrieval_param['line_file'] = line_file

retrieval_param['boxcar_size'] = 128

# known bad channel for this instrument (in the form of a numpy array):
retrieval_param['bad_channels'] = np.arange(0, 104)

retrieval_param['Tb_max'] = 200
retrieval_param['Tb_min'] = 0
retrieval_param['boxcar_thresh'] = 7

fascod_atmosphere = 'midlatitude-summer'
retrieval_param['prefix_atm'] = ARTS_DATA_PATH + \
    "/planets/Earth/Fascod/{}/{}.".format(
    fascod_atmosphere, fascod_atmosphere)

level1b_dataset, meteo_ds, global_attrs_level1b = data_GROSOM.read_level1b(
    filename)

#retrieval_param = {**global_attrs_level1b, **retrieval_param}

level1b_dataset = data_GROSOM.find_bad_channels(
    level1b_dataset, retrieval_param)

# Extracting some parameters
f0 = retrieval_param['obs_freq']
cycle = retrieval_param["integration_cycle"]
ds_freq = level1b_dataset.frequencies.values
ds_num_of_channel = len(ds_freq)
# ds_Tb = level1b_dataset.Tb[cycle].values
ds_Tb_corr = level1b_dataset.Tb_corr[cycle].values

ds_bw = max(ds_freq) - min(ds_freq)

ds_df = ds_bw/(ds_num_of_channel-1)

retrieval_param["zenith_angle"] = level1b_dataset.meanAngleAntenna.values[cycle]

retrieval_param["time"] = level1b_dataset.time[cycle].values
retrieval_param['time_start'] = level1b_dataset.firstSkyTime[cycle].values
retrieval_param['time_stop'] = level1b_dataset.lastSkyTime[cycle].values
retrieval_param["lat"] = level1b_dataset.lat[cycle].values
retrieval_param["lon"] = level1b_dataset.lon[cycle].values

retrieval_param["f_max"] = max(ds_freq)
retrieval_param["f_min"] = min(ds_freq)
retrieval_param["bandwidth"] = max(ds_freq) - min(ds_freq)

ac = arts.ArtsController()
ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)


#TODO
# Simulation grids
lat_grid = np.array([retrieval_param["lat"]])
lon_grid = np.array([retrieval_param["lon"]])

# frequency grid
n_f = retrieval_param["number_of_freq_points"]  # Number of points
bw = 1.5e9  # Bandwidth
x = np.linspace(-0.5, 0.5, n_f)
f_grid = x ** 3 + x / 10
f_grid = f_grid * bw / (max(f_grid) - min(f_grid))

f_grid = f_grid + f0

# Pressure grid != retrieval grid 
z_bottom = 1e3
z_top = 112e3
z_res = 1e3
z_grid = np.arange(z_bottom, z_top + 2 * z_res, z_res)
p_grid = z2p_simple(z_grid)

# p_grid = np.logspace(5, -1, 161)
ac.set_grids(f_grid, p_grid)

# altitude for the retrieval
retrieval_param["surface_altitude"] = 1500
ac.set_surface(retrieval_param["surface_altitude"])

# Spectroscopy
ac.set_spectroscopy_from_file(
    abs_lines_file=retrieval_param['line_file'],
    abs_species=["O3", "H2O-PWR98", "O2-PWR98", "N2-SelfContStandardType"],
    format='Arts',
    line_shape=["VVH", 750e9],
)

# Atmosphere (a priori)
# fascod_atmosphere = 'midlatitude-summer'
# prefix_atm = ARTS_DATA_PATH + "/planets/Earth/Fascod/{}/{}.".format(fascod_atmosphere,fascod_atmosphere)

retrieval_param['cira86_path'] = os.path.join(
    ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')

#apriori_data_GROSOM.plot_apriori_cira86(retrieval_param)

#fascod_atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param)

t1 = pd.to_datetime(retrieval_param['time_start'])
t2 = pd.to_datetime(retrieval_param['time_stop'])
extra_time_ecmwf = 6

ecmwf_store = '/home/eric/Documents/PhD/GROSOM/ECMWF'
cira86_path = retrieval_param['cira86_path']

atm = apriori_data_GROSOM.get_apriori_atmosphere_fascod_ecmwf_cira86(
    retrieval_param,
    ecmwf_store,
    cira86_path,
    t1,
    t2,
    extra_time_ecmwf
)

# ecmwf_atm = apriori_data_GROSOM.get_apriori_atmosphere(retrieval_param)

# Does not work now:
# The grid vector *retrieval pressure grid* is not covered by the
# corresponding atmospheric grid.
ac.set_atmosphere(atm, vmr_zeropadding=True)

# Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
ac.apply_hse(100e2, 0.1)  # value taken from GROMOS retrieval

# Sensor pos/los/time
obs = arts.Observation(
    za=40,
    aa=32,
    lat=46,
    lon=7,
    alt=15e3,
    time=12
)

ac.set_observations([obs])
ac.set_y([ds_Tb_corr])

# Defining our sensors
sensor = arts.SensorFFT(ds_freq, ds_df)
ac.set_sensor(sensor)

# doing the checks
ac.checked_calc()
y_FM, = ac.y_calc()

# Setup the retrieval grid 
# for GROMOS, 51 levels and 30 for SOMORA, we go for 32
z_bottom_ret = z_bottom
z_top_ret = 95e3
z_res_ret = 3e3
p_ret_grid = z2p_simple(np.arange(z_bottom_ret, z_top_ret, z_res_ret))

sx_O3 = covmat.covmat_diagonal_sparse(1e-6 * np.ones_like(p_ret_grid))
#sx_H2O = covmat.covmat_diagonal_sparse(1e-6 * np.ones_like(p_ret_grid))

# The different things we want to retrieve
ozone_ret = arts.AbsSpecies(
    'O3', p_ret_grid, lat_grid, lon_grid, sx_O3, unit='vmr')

# h2o_ret = arts.AbsSpecies('H2O', p_ret_grid, lat_grid, lon_grid, sx_H2O, unit='vmr')

# fshift_ret = arts.FreqShift(100e3, df=50e3)
# polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])

# increase variance for spurious spectra by factor
retrieval_param['increased_var_factor'] = 100
factor = retrieval_param['increased_var_factor']
retrieval_param['unit_var_y'] = 2e-2

y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)

y_var[(level1b_dataset.good_channels[cycle].values == 0)] = factor*retrieval_param['unit_var_y']

# y_var = utils.var_allan(ds_Tb) * np.ones_like(ds_Tb)
ac.define_retrieval([ozone_ret], y_var)
# ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)

# Let a priori be off by 0.5 ppm (testing purpose)
# vmr_offset = -0.5e-6
# ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)

# Run retrieval (parameter taken from MOPI)
# SOMORA is using 'lm': Levenberg-Marquardt (LM) method
ac.oem(
     method='lm',
     max_iter=10,
     stop_dx=0.01,
     lm_ga_settings=[100.0, 3.0, 5.0, 10.0, 1.0, 10.0],
     inversion_iterate_agenda=inversion_iterate_agenda,
)

if not ac.oem_converged:
    print("OEM did not converge.")
    print("OEM diagnostics: " + str(ac.oem_diagnostics))
    for e in ac.oem_errors:
        print("OEM error: " + e)
        continue

# Plots!
yf = ac.yf[0]
fig, axs = plt.subplots(2, sharex=True)
axs[0].plot((ds_freq - f0) / 1e6, ds_Tb_corr, label='observed')
axs[0].plot((ds_freq - f0) / 1e6, yf, label='fitted')
# axs[0].plot((ds_freq - f0) / 1e6, bl, label='baseline')
axs[0].set_ylabel('$T_B$ [K]')
axs[0].set_ylim((-10,50))
axs[0].legend()
axs[1].plot((ds_freq - f0) / 1e6, ds_Tb_corr -
            yf, label='observed - computed')
axs[1].legend()
axs[1].set_xlabel('f - {:.3f} GHz [MHz]'.format(f0/1e9))
axs[1].set_ylabel('$T_B$ [K]')
axs[1].set_ylim((-10,10))
fig.tight_layout()
if save_plots:
    fig.savefig(level2_data_folder+name+'_spectrum.pdf')
if show_plots:
    fig.show()
fig, axs = plt.subplots(2, 2, sharey=True)
axs[0][0].plot(ozone_ret.x * 1e6, ozone_ret.z_grid /
               1e3, label='retrieved', marker='x')
axs[0][0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label='apriori')
# axs[0][0].plot((ozone_ret.xa - vmr_offset) * 1e6, ozone_ret.z_grid / 1e3, label='true')
axs[0][0].set_xlabel('Ozone VMR [ppm]')
axs[0][0].set_ylabel('Altitude [km]')
axs[0][0].legend()
axs[0][1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
axs[0][1].set_xlabel('Measurement response')
axs[1][0].plot(ozone_ret.offset / 1e3, ozone_ret.z_grid / 1e3)
axs[1][0].set_xlabel('Offset [km]')
axs[1][0].set_ylabel('Altitude [km]')
for avk in ozone_ret.avkm:
    if 0.8 <= np.sum(avk) <= 1.2:
        axs[1][1].plot(avk, ozone_ret.z_grid / 1e3)
fig.tight_layout()

if save_plots:
    fig.savefig(level2_data_folder+name+'.pdf')
    print('\nSaved plots to :', level2_data_folder)

if show_plots:
    fig.show()

if save_netcdf:
    ac.get_level2_xarray().to_netcdf(level2_data_folder+name+'.nc')
    print('\nSaved results to :', level2_data_folder)

