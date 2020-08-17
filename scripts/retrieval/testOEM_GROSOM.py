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
from retrievals import utils
from retrievals.arts.atmosphere import p2z_simple, z2p_simple
from retrievals.data import p_interpolate
#from pyarts.workspace import arts_agenda
from typhon.arts.workspace import arts_agenda
from dotenv import load_dotenv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings('ignore', message='numpy.dtype size changed')

load_dotenv('/home/esauvageat/Documents/ARTS/arts-examples/.env')
#load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

show_plots = True
save_plots = True
save_netcdf = False

name = 'testOEM_SOMORA_daily_Fascod'

# For testing
#basename = "/home/eric/Documents/PhD/GROSOM/Level1/"
#level2_data_folder = "/home/eric/Documents/PhD/GROSOM/Level2/"

basename="/scratch/GROSOM/Level1/SOMORA/"
level2_data_folder = "/scratch/GROSOM/Level2/"

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
    ws.atmfields_checkedCalc(negative_vmr_ok=1)
    ws.atmgeom_checkedCalc()

    # Calculate yf and Jacobian matching x
    ws.yCalc(y=ws.yf)

    # Add baseline term
    ws.VectorAddVector(ws.yf, ws.y, ws.y_baseline)

    # This method takes cares of some "fixes" that are needed to get the Jacobian
    # right for iterative solutions. No need to call this WSM for linear inversions.
    ws.jacobianAdjustAndTransform()

instrument_name = "SOMORA"
filename = basename+"SOMORA_level1b_24h_AC240_2019_02_04"

retrieval_param = dict()
retrieval_param["integration_cycle"] = 0
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

level1b_dataset, flags, meteo_ds, global_attrs_level1b = data_GROSOM.read_level1b(
    filename)

#retrieval_param = {**global_attrs_level1b, **retrieval_param}
cycle = retrieval_param["integration_cycle"]
#level1b_dataset = data_GROSOM.find_bad_channels(
#    level1b_dataset, np.arange(0, 104), 0, 260, 128, 7)

#goodFreq = np.where(np.logical_and(level1b_dataset.intermediate_freq.data <=1000,level1b_dataset.stdTb[cycle].data < 20))[0]
goodFreq = np.where(level1b_dataset.intermediate_freq.data <=1000)
# Extracting some parameters
f0 = retrieval_param['obs_freq']

ds_freq = level1b_dataset.frequencies.values[goodFreq]
ds_num_of_channel = len(ds_freq)
ds_Tb = level1b_dataset.Tb[cycle].values
ds_Tb_corr = level1b_dataset.Tb_corr[cycle].values
ds_Tb_corr_clean = ds_Tb_corr[goodFreq]

ds_bw = max(ds_freq) - min(ds_freq)

ds_df = ds_bw/(ds_num_of_channel-1)

retrieval_param["zenith_angle"] = 90 - level1b_dataset.mean_sky_elevation_angle.values[cycle]
retrieval_param["azimuth_angle"] = level1b_dataset.azimuth_angle.values[cycle]
retrieval_param["time"] = level1b_dataset.time[cycle].values
retrieval_param["lat"] = level1b_dataset.lat[cycle].values
retrieval_param["lon"] = level1b_dataset.lon[cycle].values
retrieval_param['time_start'] = level1b_dataset.first_sky_time[cycle].values
retrieval_param['time_stop'] = level1b_dataset.last_sky_time[cycle].values
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

bw = 1.3e9  # Bandwidth
x = np.linspace(-1, 1, n_f)
f_grid = x ** 3 + x / 50
f_grid = f_grid * bw / (max(f_grid) - min(f_grid))
f_grid = f_grid + f0
# Trying with another frequency grid
#centerIF = [475e6, 525e6]
#centerGrid = np.arange(retrieval_param["f_min"] + centerIF[0], retrieval_param["f_min"] + centerIF[1],ds_df+1e3)

#leftGrid = np.arange(retrieval_param["f_min"], retrieval_param["f_min"] + centerIF[0], 1e6)
#rightGrid = np.arange(retrieval_param["f_min"]+centerIF[1], retrieval_param["f_max"]+1e6, 1e6)

#f_grid = np.hstack((leftGrid, centerGrid,rightGrid))


# Pressure grid != retrieval grid 
z_bottom = 800
z_top = 95e3
z_res = 1e3
z_grid = np.arange(z_bottom, z_top + 2 * z_res, z_res)
p_grid = z2p_simple(z_grid)

#p_grid = np.logspace(5, -1, 161)
ac.set_grids(f_grid, p_grid)

# altitude for the retrieval
retrieval_param["surface_altitude"] = 10e3
ac.set_surface(retrieval_param["surface_altitude"])

# Spectroscopy
ac.set_spectroscopy_from_file(
    abs_lines_file=retrieval_param['line_file'],
    abs_species=["O3", "H2O-PWR98", "O2-PWR93", "N2-SelfContStandardType"],
    format='Arts',
    line_shape=["VVH", 750e9],
)

# Atmosphere (a priori)
# fascod_atmosphere = 'midlatitude-summer'
# prefix_atm = ARTS_DATA_PATH + "/planets/Earth/Fascod/{}/{}.".format(fascod_atmosphere,fascod_atmosphere)
"""
retrieval_param['cira86_path'] = os.path.join(
    ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')

#apriori_data_GROSOM.plot_apriori_cira86(retrieval_param)

#fascod_atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param)

t1 = pd.to_datetime(retrieval_param['time_start'])
t2 = pd.to_datetime(retrieval_param['time_stop'])
extra_time_ecmwf = 4

#ecmwf_store = '/home/eric/Documents/PhD/GROSOM/ECMWF'
ecmwf_store = '/home/esauvageat/Documents/GROSOM/Analysis/ECMWF'
cira86_path = retrieval_param['cira86_path']

atm = apriori_data_GROSOM.get_apriori_atmosphere_fascod_ecmwf_cira86(
    retrieval_param,
    ecmwf_store,
    cira86_path,
    t1,
    t2,
    extra_time_ecmwf
)

#ecmwf_atm = apriori_data_GROSOM.get_apriori_atmosphere(retrieval_param)

# Does not work now:
# The grid vector *retrieval pressure grid* is not covered by the
# corresponding atmospheric grid.
ac.set_atmosphere(atm, vmr_zeropadding=True)

"""

ac.set_atmosphere_fascod('midlatitude-winter')

# Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
ac.apply_hse(100e2, 0.5)  # value taken from MOPI retrieval

# Sensor pos/los/time
obs = arts.Observation(
    za = retrieval_param["zenith_angle"], 
    aa = retrieval_param["azimuth_angle"], 
    lat = retrieval_param["lat"],
    lon = retrieval_param["lon"],
    alt = retrieval_param["observation_altitude"],
    time = retrieval_param["time"]
    )

ac.set_observations([obs])
ac.set_y([ds_Tb_corr_clean])

# Defining our sensors
sensor = arts.SensorFFT(ds_freq, ds_df)
ac.set_sensor(sensor)

# doing the checks
ac.checked_calc()
#y_FM, = ac.y_calc()

# Setup the retrieval grid 
# for GROMOS, 51 levels and 30 for SOMORA, we go for 32
z_bottom_ret = z_bottom
z_top_ret = 95e3
z_res_ret = 3e3
p_ret_grid = z2p_simple(np.arange(z_bottom_ret, z_top_ret, z_res_ret))
#p_ret_grid = np.logspace(5, -1, 32)

'''
o3_std_value=0.7e-6
o3_std_const = o3_std_value * np.ones_like(p_ret_grid)
o3_apriori = atm.vmr_field('o3').data[:,0][:,0]
o3_std_limit = p_interpolate(
    p_ret_grid, atm.vmr_field('o3').grids[0], o3_apriori, fill=0
)

o3_std = np.minimum(o3_std_limit, o3_std_const)
#std_correction_o3 = 1
ozone_covmat = covmat.covmat_1d_sparse(
    np.log10(p_ret_grid),
    1e-6 * np.ones_like(p_ret_grid),
    0.3 * np.ones_like(p_ret_grid),
    fname="lin",
    cutoff=0.001,
)

h2o_std_value=0.7e-6
h2o_std_const = h2o_std_value * np.ones_like(p_ret_grid)
h2o_apriori = atm.vmr_field('h2o').data[:,0][:,0]
h2o_std_limit = p_interpolate(
    p_ret_grid, atm.vmr_field('h2o').grids[0], h2o_apriori, fill=0
)

h2o_std = np.minimum(h2o_std_limit, h2o_std_const)
std_correction_h2o = 1
h2o_covmat = covmat.covmat_1d_sparse(
    np.log10(p_ret_grid),
    std_correction_h2o * h2o_std,
    0.3 * np.ones_like(p_ret_grid),
    fname="lin",
    cutoff=0.001,
)
'''
sx_O3 = covmat.covmat_diagonal_sparse(1.5e-6 * np.ones_like(p_ret_grid))
#sx_H2O = covmat.covmat_diagonal_sparse(1e-12 * np.ones_like(p_ret_grid))

fshift_ret = arts.FreqShift(100e3, df=50e3)

# Load the a priori atmospheric state.
#from typhon.arts import xml
#atm_fields = xml.load("/home/esauvageat/Documents/ARTS/arts-lectures/exercises/06-inversion/input/x_apriori.xml")
#z = atm_fields.get("z", keep_dims=False)
#x_apriori = atm_fields.get("abs_species-H2O", keep_dims=False)

# Load the covariance matrices.
#sx_H2O_artsLectures = xml.load("/home/esauvageat/Documents/ARTS/arts-lectures/exercises/06-inversion/input/S_xa.xml")

# The different things we want to retrieve
ozone_ret = arts.AbsSpecies('O3', p_ret_grid, lat_grid, lon_grid, sx_O3, unit='vmr')
#h2o_ret = arts.AbsSpecies("H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252", p_ret_grid, lat_grid, lon_grid, sx_H2O, unit='vmr')

# fshift_ret = arts.FreqShift(100e3, df=50e3)
# polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])

# increase variance for spurious spectra by factor
#retrieval_param['increased_var_factor'] = 100
#factor = retrieval_param['increased_var_factor']
#retrieval_param['unit_var_y'] = 0.1

#y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_Tb_corr_clean)
# var = np.mean(np.square(np.diff(ds_Tb_corr))) / 2
#y_var = 10*utils.var_allan(ds_Tb_corr) * np.ones_like(ds_Tb_corr)

#var = np.mean(np.square(np.diff(ds_Tb_corr, axis=0)), axis=0) / 2
#var = np.mean(np.square(np.diff(ds_Tb_corr_clean, axis=0)), axis=0) / 2

#y_var = 100*var * np.ones_like(ds_Tb_corr_clean)

y_var = 10*level1b_dataset.stdTb[cycle].data[goodFreq]
#y_var[(level1b_dataset.good_channels[cycle].values == 0)] = factor*retrieval_param['unit_var_y']
#polyfit_ret = arts.Polyfit(
#    poly_order=1, covmats=[np.array([[5]]), np.array([[1]])]
#)

ac.define_retrieval([ozone_ret, fshift_ret], y_var)
# ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)

#Let a priori be off by 0.5 ppm (testing purpose)
#vmr_offset = +0.5e-6
#ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)

# Run retrieval (parameter taken from MOPI)
# SOMORA is using 'lm': Levenberg-Marquardt (LM) method
# Run retrieval
ac.oem(
    method="gn",
    max_iter=10,
    stop_dx=0.1,
    inversion_iterate_agenda=inversion_iterate_agenda,
)
    
if not ac.oem_converged:
    print("OEM did not converge.")
    print("OEM diagnostics: " + str(ac.oem_diagnostics))
    for e in ac.oem_errors:
        print("OEM error: " + e)
        continue

figures = list()

# Plots!
yf = ac.yf[0]
r = ds_Tb_corr_clean - yf
r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")

fig, axs = plt.subplots(2, sharex=True)
axs[0].plot((ds_freq - f0) / 1e6, ds_Tb_corr_clean, label='observed')
axs[0].plot((ds_freq - f0) / 1e6, yf, label='fitted')
#axs[0].plot((f_grid - f0) / 1e6, 40*np.ones(len(f_grid)), '.', label='grid points')
axs[0].set_ylabel('$T_B$ [K]')
axs[0].set_ylim((-10,50))
#axs[0].set_xlim((-25,25))
axs[0].legend()
axs[1].plot((ds_freq - f0) / 1e6, r, label='observed - computed')
axs[1].plot((ds_freq - f0) / 1e6, r_smooth, label="residuals smooth")
axs[1].legend()
axs[1].set_xlabel('f - {:.3f} GHz [MHz]'.format(f0/1e9))
axs[1].set_ylabel('$T_B$ [K]')
axs[1].set_ylim((-2,2))
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
figures.append(fig)
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
axs[1][0].plot(ozone_ret.es * 1e6, ozone_ret.z_grid / 1e3, label="smoothing error")
axs[1][0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3, label="obs error")
#axs[1][0].plot(ozone_ret.offset / 1e3, ozone_ret.z_grid / 1e3)
axs[1][0].set_xlabel('Error [ppm]')
axs[1][0].set_ylabel('Altitude [km]')
for avk in ozone_ret.avkm:
    if 0.8 <= np.sum(avk) <= 1.2:
        axs[1][1].plot(avk, ozone_ret.z_grid / 1e3)

axs[0][0].grid(True)
axs[0][1].grid(True)
axs[1][1].grid(True)
axs[1][0].grid(True)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
figures.append(fig)

if save_plots:
    filename = level2_data_folder+name+'.pdf'
    with PdfPages(filename) as pdf:
        for fig in figures:
            pdf.savefig(fig)
    #fig.savefig(level2_data_folder+name+'.pdf')
    print('\nSaved plots as :', filename)

if show_plots:
    fig.show()

if save_netcdf:
    ac.get_level2_xarray().to_netcdf(level2_data_folder+name+'.nc')
    print('\nSaved results to :', level2_data_folder)