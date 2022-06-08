"""
This is very close to the "controlfiles/artscomponents/TestOEM.arts" cfile shipped with arts.
Plus some plots to show the retrieval results.
"""
import sys
sys.path.insert(0, '/home/es19m597/Documents/GROMORA/GROMORA-harmo/scripts/pyretrievals/')
import warnings
warnings.filterwarnings('ignore', message='numpy.dtype size changed')
import os
import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv

from retrievals import arts
from retrievals import covmat

from pyarts.workspace import arts_agenda

# For ARTS, we need to specify some paths
load_dotenv('/opt/anaconda/.env.birg-arts24')
#load_dotenv('/opt/arts/.env.stockhorn-arts24')

ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

def make_f_grid():
    n_f = 601  # Number of points
    bw = 600e6  # Bandwidth
    x = np.linspace(-1, 1, n_f)
    f_grid = x ** 3 + x / 10
    f_grid = f_grid * bw / (max(f_grid) - min(f_grid))
    return f_grid

@arts_agenda
def gromora_inversion_agenda(ws):
    """Custom inversion iterate agenda to ignore bad partition functions."""
    ws.Ignore(ws.inversion_iteration_counter)

    ws.xClip(ijq=0, limit_low=0.00000000001, limit_high=0.00002)

    # Map x to ARTS' variables
    ws.x2artsAtmAndSurf()
    ws.x2artsSensor()

    # To be safe, rerun some checkss
    ws.atmfields_checkedCalc(negative_vmr_ok=True)
    ws.atmgeom_checkedCalc()

    # Calculate yf and Jacobian matching x
    ws.yCalc() #()ws.yf

    # Add baseline term
    #ws.VectorAddElementwise(ws.yf, ws.y, ws.y_baseline)
    ws.VectorAddVector(ws.yf, ws.y, ws.y_baseline)

    # This method takes cares of some "fixes" that are needed to get the Jacobian
    # right for iterative solutions. No need to call this WSM for linear inversions.
    ws.jacobianAdjustAndTransform()

def main(show_plots=False, save_plots=False, save_netcdf=False):
    f0 = 110.836e+9

    ac = arts.ArtsController()
    ac.setup(atmosphere_dim=1)

    # Simulation grids
    f_grid = make_f_grid() + f0
    p_grid = np.logspace(5, -1, 361)
    ac.set_grids(f_grid, p_grid)

    # Spectroscopy
    ac.set_spectroscopy_from_file(
            abs_lines_file='spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz',
            abs_species=['O3'],
            format='Arts',
            line_shape=["VVH", 750e9],
        )
    # ac.set_spectroscopy_from_file('spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz', ['O3'],
    #                               'Arts', ('Voigt_Kuntz6', 'VVH', 750e9))

    # Atmosphere (a priori)
    clim_prefix = os.path.join(ARTS_DATA_PATH, 'planets/Earth/Fascod/tropical/tropical.')
    atmosphere = arts.Atmosphere.from_arts_xml(clim_prefix)
    ac.set_atmosphere(atmosphere)
    ac.set_surface(10e3)

    # Apply HSE
    ac.apply_hse(100e2, 0.5)

    # Sensor pos/los/time
    obs1 = arts.Observation(za=60, aa=0, alt=15e3)
    ac.set_observations([obs1, ])

    # Backend
    f_resolution = 200e3
    f_backend = np.arange(-280e6, 280e6, f_resolution) + f0
    sensor = arts.SensorGaussian(f_backend, np.array([f_resolution]))
    ac.set_sensor(sensor)

    # Perform checks
    ac.checked_calc()
    y, = ac.y_calc()  # only one observation

    # Setup retrieval
    p_ret_grid = np.logspace(5, -1, 81)
    lat_ret_grid = lon_ret_grid = np.array([0])
    sx = covmat.covmat_diagonal_sparse(1e-6 * np.ones_like(p_ret_grid))
    ozone_ret = arts.AbsSpecies('O3', p_ret_grid, lat_ret_grid, lon_ret_grid, sx, unit='vmr')

    fshift_ret = arts.FreqShift(100e3, df=50e3)

    polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])

    y_var = 1e-2 * np.ones_like(f_backend)

    ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)

    # Let a priori be off by 0.5 ppm
    vmr_offset = 0.5e-6
    ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)

    # Add a baseline
    bl_coeffs = [2, 1]
    y = y + bl_coeffs[0] * np.linspace(0, 1, y.size) + bl_coeffs[1]
    ac.set_y([y, ])

    # Use shifted sensor
    f_shift = -150e3
    ac.set_sensor(arts.SensorGaussian(f_backend + f_shift, np.array([f_resolution])))

    # Run retrieval
    ac.oem(method='gn', max_iter=10, stop_dx=0.1, inversion_iterate_agenda=gromora_inversion_agenda)

    if not ac.oem_converged:
        return False

    # Plots!
    yf = ac.yf[0]
    bl = ac.y_baseline[0]
    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - f0) / 1e6, y, label='observed')
    axs[0].plot((f_backend - f0) / 1e6, yf, label='fitted')
    axs[0].plot((f_backend - f0) / 1e6, bl, label='baseline')
    axs[0].legend()
    axs[1].plot((f_backend - f0) / 1e6, y-yf, label='observed - computed')
    axs[1].legend()
    axs[1].set_xlabel('f - {:.3f} GHz [MHz]'.format(f0/1e9))
    axs[0].set_ylabel('$T_B$ [K]')
    axs[1].set_ylabel('$T_B$ [K]')
    fig.tight_layout()
    if save_plots:
        fig.savefig('test_oem_spectrum.png')
    if show_plots:
        fig.show()

    fig, axs = plt.subplots(2, 2, sharey=True)

    axs[0][0].plot(ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label='retrieved', marker='x')
    axs[0][0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label='apriori')
    axs[0][0].plot((ozone_ret.xa - vmr_offset) * 1e6, ozone_ret.z_grid / 1e3, label='true')
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
        fig.savefig('test_oem_ozone.png')
        print('\nSaved plots to TestOEM_*.png')
    if show_plots:
        fig.show()

    print('Fshift fit: {:g} kHz, true: {:g} kHz'.format(fshift_ret.x[0]/1e3, f_shift/1e3))
    print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x])
          + ' true: '+ ', '.join(map(str, bl_coeffs)))

    if save_netcdf:
        ac.get_level2_xarray().to_netcdf('TestOEM_result.nc')
        print('\nSaved results to TestOEM_result.nc')
    
    return True


def test_oem():
    """ Run this example as test """
    assert main(show_plots=True, save_plots=False, save_netcdf=False)


if __name__ == '__main__':
    main(show_plots=True, save_plots=False, save_netcdf=False)

