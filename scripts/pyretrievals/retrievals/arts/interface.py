import numpy as np
import xarray as xr
import scipy
from scipy import sparse
from typing import NamedTuple
import os
import uuid
import datetime
#from dotenv import load_dotenv

#load_dotenv(dotenv_path='./.env')
import matplotlib.pyplot as plt
#from typhon.arts.workspace import Workspace, arts_agenda
#from typhon.arts.workspace import Workspace
from pyarts.workspace import Workspace
#from pyarts.workspace import arts_agenda
from pyarts import xml
from pyarts import classes
from pyarts.griddedfield import GriddedField3

from retrievals.arts import boilerplate
#from retrievals.arts import retrieval
from retrievals.data import p_interpolate
from retrievals.data.utils import interpolate


def _is_asc(x):
    """Check if vector is strongly monotonic increasing."""
    return np.all(np.diff(x) > 0)


def _is_desc(x):
    """Check if vector is strongly monotonic decreasing."""
    return np.all(np.diff(x) < 0)


class Observation(NamedTuple):
    """
    Geometry related to an observation.
    """

    #: Zenith angle in [0, 180] (for 1D and 3D atmosphere). See :arts:variable:`sensor_los`.
    za: float

    #: Azimuth angle in [-180, 180] counted clockwise, 0 is north and 90 is east. See :arts:variable:`sensor_los`.
    aa: float

    #: Latitude of sensor, default 0. See :arts:variable:`sensor_pos`.
    lat: float = 0

    #: Longitude of sensor, default 0. See :arts:variable:`sensor_pos`.
    lon: float = 0

    #: Altitude of sensor, default 0. See :arts:variable:`sensor_pos`.
    alt: float = 0

    #: Time of sensor, default 0. See :arts:variable:`sensor_time`.
    time: float = 0


class OemException(Exception):
    """Wrapper for all exceptions that happen during a OEM WSM call."""
    pass


class ArtsController():
    """A not so high level interface to ARTS."""

    def __init__(self, verbosity=0, agenda_verbosity=0):
        self.ws = Workspace(verbosity=verbosity, agenda_verbosity=agenda_verbosity)
        self.retrieval_quantities = []
        self._sensor = None
        self._observations = []

    def setup(self, atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1):
        """
        Run boilerplate (includes, agendas) and set basic variables.

        :param atmosphere_dim:
        :param iy_unit:
        :param ppath_lmax:
        :param stokes_dim:
        """
        boilerplate.include_general(self.ws)
        boilerplate.copy_agendas(self.ws)
        boilerplate.set_basics(self.ws, atmosphere_dim, iy_unit, ppath_lmax, stokes_dim)

        # Deactivate some stuff (can be activated later)
        self.ws.jacobianOff()
        self.ws.cloudboxOff()

    def checked_calc(self, negative_vmr_ok=False):
        """Run checked calculations."""
        boilerplate.run_checks(self.ws, negative_vmr_ok)

    def set_spectroscopy(self, abs_lines, abs_species, line_shape=None, abs_f_interp_order=3):
        """
        Setup absorption species and spectroscopy data.

        :param ws: The workspace.
        :param abs_lines: Absoption lines.
        :param abs_species: List of abs species tags.
        :param line_shape: Line shape definition. Default: ['Voigt_Kuntz6', 'VVH', 750e9]
        :param f_abs_interp_order: No effect for OnTheFly propmat. Default: 3
        :type abs_lines: typhon.arts.catalogues.ArrayOfLineRecord
        """
        boilerplate.setup_spectroscopy(self.ws, abs_lines, abs_species, line_shape)
        self.ws.abs_f_interp_order = abs_f_interp_order  # no effect for OnTheFly propmat
    
    def set_spectroscopy_from_file_old(self, abs_lines_file, abs_species,  format='HITRAN', line_shape=None, abs_f_interp_order=3):
        """
        Setup absorption species and spectroscopy data from HITRAN catalogue file.

        :param ws: The workspace.
        :param basename:  Path to an XML file.
        :param abs_species: List of abs species tags.
        :param format: One of 'ARTSCAT', 'JPL', 'HITRAN' (and others for which a WSM `Read...` exists)
        :param line_shape: Line shape definition. Default: ['VVH', 750e9]
        :param f_abs_interp_order: No effect for OnTheFly propmat. Default: 3
        """
        ws = self.ws
        if line_shape is None:
            line_shape = ['VVH', 750e9]
        
        ws.abs_speciesSet(abs_species)
        #ws.abs_lineshapeDefine(*line_shape)
        
        #TODO make it work for all types
        read_fn = getattr(ws, 'Read' + format)
        read_fn(ws.abs_lines, filename = abs_lines_file)
        
        ws.abs_linesSetNormalization(ws.abs_lines, line_shape[0])
        ws.abs_linesSetCutoff(ws.abs_lines, "ByLine", line_shape[1])

        ws.abs_lines_per_speciesCreateFromLines()
        ws.abs_f_interp_order = abs_f_interp_order    
    
    def set_spectroscopy_from_file(self, abs_lines_file, abs_species, format='Arts', line_shape=None, abs_f_interp_order=3):
        """
        Setup absorption species and spectroscopy data from XML file.

        :param ws: The workspace.
        :param abs_lines_file: Path to an XML file.
        :param abs_species: List of abs species tags.
        :param format: One of 'Arts', 'Jpl', 'Hitran' (and others for which a WSM `abs_linesReadFrom...` exists)
        :param line_shape: Line shape definition. Default: ['VVH', 750e9]
        :param f_abs_interp_order: No effect for OnTheFly propmat. Default: 3
        """
        ws = self.ws
        if line_shape is None:
            line_shape = ['VVH', 750e9]
        ws.abs_speciesSet(species=abs_species)
        #ws.abs_lineshapeDefine(*line_shape)
        
        ws.ReadXML(ws.abs_lines, abs_lines_file)
        
        ws.abs_linesSetCutoff(ws.abs_lines, "ByLine", line_shape[1])
        ws.abs_linesSetNormalization(ws.abs_lines, line_shape[0])
        
        ws.abs_lines_per_speciesCreateFromLines()
        ws.abs_f_interp_order = abs_f_interp_order

    def set_grids(self, f_grid, p_grid, lat_grid=None, lon_grid=None):
        """
        Set the forward model grids. Basic checks are performed, depending on dimensionality of atmosphere.

        :param f_grid:
        :param p_grid:
        :param lat_grid:
        :param lon_grid:
        """
        if lat_grid is None:
            lat_grid = np.array([0])
        if lon_grid is None:
            lon_grid = np.array([0])

        if not self.ws.atmosphere_dim.initialized:
            raise Exception('atmosphere_dim must be initialized before assigning grids.')
        if not _is_asc(f_grid):
            raise ValueError('Values of f_grid must be strictly increasing.')
        if not _is_desc(p_grid) or not np.all(p_grid > 0):
            raise ValueError('Values of p_grid must be strictly decreasing and positive.')
        if self.ws.atmosphere_dim.value == 1:
            if lat_grid.size > 1 or lon_grid.size > 1:
                raise ValueError('For 1D atmosphere, lat_grid and lon_grid shall be of length 1.')
        elif self.ws.atmosphere_dim.value == 2:
            if lon_grid is not None or len(lat_grid):
                raise ValueError('For 2D atmosphere, lon_grid shall be empty.')
            if lat_grid is None or len(lat_grid) == 0:
                raise ValueError('For 2D atmosphere, lat_grid must be set.')
        elif self.ws.atmosphere_dim.value == 3:
            if lat_grid is None or len(lat_grid) < 2 or lon_grid is None or len(lon_grid) < 2:
                raise ValueError('For 3D atmosphere, length of lat_grid and lon_grid must be >= 2.')
            if max(abs(lon_grid)) > 360:
                raise ValueError('Values of lon_grid must be in the range [-360,360].')
            if max(abs(lat_grid)) > 90:
                raise ValueError('Values of lat_grid must be in the range [-90,90].')
        if lat_grid is not None and not _is_asc(lat_grid):
            raise ValueError('Values of lat_grid must be strictly increasing.')
        if lon_grid is not None and not _is_asc(lon_grid):
            raise ValueError('Values of lon_grid must be strictly increasing.')

        self.ws.f_grid = f_grid
        self.ws.p_grid = p_grid
        if self.atmosphere_dim > 1:
            self.ws.lat_grid = lat_grid
            self.ws.lon_grid = lon_grid
        else:
            self.ws.lat_grid = []
            self.ws.lon_grid = []
            self.ws.lat_true = lat_grid
            self.ws.lon_true = lon_grid
        self.set_surface(0)

    def set_surface(self, altitude: float):
        """
        Set surface altitude.

        :param float altitude:

        .. note:: Currently, only constant altitudes are supported.
        """
        shape = max((1, 1), (self.n_lat, self.n_lat))
        self.ws.z_surface = altitude * np.ones(shape)

    def set_atmosphere(self, atmosphere, vmr_zeropadding=False):
        """
        Set the atmospheric state.

        :param atmosphere: Atmosphere with Temperature, Altitude and VMR Fields.
        :type atmosphere: retrievals.arts.atmosphere.Atmosphere
        :param vmr_zeropadding: Allow VMR zero padding wind fields are always zero
                                padded. Default: False.

        .. note:: Currently only supports 1D atmospheres that is then expanded to a
            multi-dimensional homogeneous atmosphere.
        """
        vmr_zeropadding = 1 if vmr_zeropadding else 0

        self.ws.t_field_raw = atmosphere.t_field
        self.ws.z_field_raw = atmosphere.z_field
        self.ws.vmr_field_raw = [atmosphere.vmr_field(mt) for mt in self.abs_species_maintags]
        self.ws.nlte_field_raw = None
        # fix for arts2.3
        self.ws.nlte_vibrational_energies = []

        for c in ('u', 'v', 'w'):
            try:
                raw_field = atmosphere.wind_field(c)
                field = p_interpolate(
                    self.p_grid, raw_field.grids[0], raw_field.data[:, 0, 0], fill=0
                )
            except KeyError:
                field = np.zeros_like(self.p_grid)
            field = np.tile(
                field[:, np.newaxis, np.newaxis], (1, self.n_lat, self.n_lon)
            )
            field_name = 'wind_{}_field'.format(c)
            setattr(self.ws, field_name, field)

        if  self.atmosphere_dim == 1:
            self.ws.AtmFieldsCalc(vmr_zeropadding=vmr_zeropadding)
        else:
            self.ws.AtmFieldsCalcExpand1D(vmr_zeropadding=vmr_zeropadding)

    def apply_hse(self, p_hse=100e2, z_hse_accuracy=0.5):
        """
        Calculate z field from hydrostatic equilibrium. See :arts:method:`z_fieldFromHSE`.

        :param p_hse: See :arts:variable:`p_hse`.
        :param z_hse_accuracy: See :arts:variable:`z_hse_accuracy`.
        """
        ws = self.ws
        ws.p_hse = p_hse
        ws.z_hse_accuracy = z_hse_accuracy
        ws.atmfields_checkedCalc()
        ws.z_fieldFromHSE()

    def set_wind(self, wind_u=None, wind_v=None, wind_w=None):
        """
        Set the wind fields to constant values.
        """
        ws = self.ws
        n_p, n_lat, n_lon = self.n_p, self.n_lat, self.n_lon

        if wind_u is not None:
            ws.Tensor3SetConstant(ws.wind_u_field, n_p, n_lat, n_lon, float(wind_u))
        if wind_v is not None:
            ws.Tensor3SetConstant(ws.wind_v_field, n_p, n_lat, n_lon, float(wind_v))
        if wind_w is not None:
            ws.Tensor3SetConstant(ws.wind_w_field, n_p, n_lat, n_lon, float(wind_w))

    def _check_set_wind(self):
        """Set uninitialized wind fields to 0."""
        ws = self.ws
        if ws.wind_u_field.value.size == 0:
            self.set_wind(wind_u=0)
        if ws.wind_v_field.value.size == 0:
            self.set_wind(wind_v=0)
        if ws.wind_w_field.value.size == 0:
            self.set_wind(wind_w=0)

    def set_observations(self, observations):
        """
        Set the geometry of the observations made.

        :param observations:
        :type observations: Iterable[retrievals.arts.interface.Observation]
        """
        sensor_los = np.array([[obs.za, obs.aa] for obs in observations])
        sensor_pos = np.array([[obs.alt, obs.lat, obs.lon] for obs in observations])
        self.ws.sensor_los = sensor_los[:, :max(self.atmosphere_dim - 1, 1)]
        self.ws.sensor_pos = sensor_pos[:, :self.atmosphere_dim]
        #self.ws.sensor_time = np.array([1 for obs in observations])#for obs in observations
        classes.from_workspace(self.ws.sensor_time).data = np.array([obs.time for obs in observations])
        #self.ws.sensor_time = np.array([obs.time for obs in observations])
        self._observations = observations

    def set_y(self, ys):
        """
        Set the observations.

        :param ys: List with a spectrum for every observation.
        :return:
        """
        y = np.concatenate(ys)
        self.ws.y = y

    def set_sensor(self, sensor):
        """
        Set the sensor.

        :param sensor:
        :type sensor: retrievals.arts.sensors.AbstractSensor
        """
        self._sensor = sensor
        sensor.apply(self.ws)

    def y_calc(self, jacobian_do=False):
        """
        Run the forward model.

        :param jacobian_do: Not implemented yet.
        :return: The measurements as list with length according to observations.
        """
        if jacobian_do:
            raise NotImplementedError('Jacobian not implemented yet.')
            # self.ws.jacobian_do = 1
        self.ws.yCalc()
        return self.y

    def define_retrieval(self, retrieval_quantities, y_vars):
        """
        Define the retrieval quantities.

        :param retrieval_quantities: Iterable of retrieval quantities `retrievals.arts.retrieval.RetrievalQuantity`.
        :param y_vars: List or variance vectors according.
        """
        ws = self.ws

        if isinstance(y_vars, (list, tuple)):
            y_vars = np.concatenate(y_vars)

        if len(y_vars) != self.n_y:
            raise ValueError('Variance vector y_vars must have same length as y.')
        
        ws.retrievalDefInit()

        # Retrieval quantities
        self.retrieval_quantities = retrieval_quantities
        for rq in retrieval_quantities:
            rq.apply(ws)

        # Se and its inverse
        covmat_block = sparse.diags(y_vars, format='csr')
        boilerplate.set_variable_by_xml(ws, ws.covmat_block, covmat_block)
        ws.covmat_seAddBlock(block=ws.covmat_block)

        covmat_block = sparse.diags(1/y_vars, format='csr')
        boilerplate.set_variable_by_xml(ws, ws.covmat_block, covmat_block)
        ws.covmat_seAddInverseBlock(block=ws.covmat_block)

        ws.retrievalDefClose()

    def oem(self, method='li', max_iter=10, stop_dx=0.01, lm_ga_settings=None, display_progress=True,
            inversion_iterate_agenda=None):
        """
        Run the optimal estimation. See Arts documentation for details.

        :param method:
        :param max_iter:
        :param stop_dx:
        :param lm_ga_settings: Default: [10, 2, 2, 100, 1, 99]
        :param display_progress:
        :param inversion_iterate_agenda: If set to None, a simple default agenda is used.
        :return:
        """
        ws = self.ws

        if lm_ga_settings is None:
            lm_ga_settings = [100.0, 2.0, 2.0, 10.0, 1.0, 1.0]
        lm_ga_settings = np.array(lm_ga_settings)

        # x, jacobian and yf must be initialised
        ws.x = np.array([])
        ws.yf = np.array([])
        ws.jacobian = np.array([[]])

        if inversion_iterate_agenda is None:
            inversion_iterate_agenda = boilerplate.inversion_iterate_agenda
        ws.Copy(ws.inversion_iterate_agenda, inversion_iterate_agenda)

        ws.AgendaExecute(ws.sensor_response_agenda)

        # a priori values
        self._check_set_wind()
        ws.xaStandard()
        xa = ws.xa.value
        for rq in self.retrieval_quantities:
            rq.extract_apriori(xa)

        #ws.covmat_sxExtractSqrtDiagonal(ws, ws.covmat_sx)
        #x_norm = np.sqrt(np.diag(ws.covmat_sx.value.to_dense()))

        # Run inversion
        try:
            ws.OEM(method=method, max_iter=max_iter, stop_dx=stop_dx, lm_ga_settings=lm_ga_settings,
                   display_progress=1 if display_progress else 0)
            # ws.OEM(method=method, x_norm=x_norm, max_iter=max_iter, stop_dx=stop_dx, lm_ga_settings=lm_ga_settings,
            #        display_progress=1 if display_progress else 0)
        except Exception as e:
            raise OemException(e)

        if not self.oem_converged:  # Just checks if dxdy is initialized
            return False

        ws.x2artsAtmAndSurf()
        ws.x2artsSensor()
        ws.avkCalc()
        ws.covmat_ssCalc()
        ws.covmat_soCalc()
        ws.retrievalErrorsExtract()

        x = ws.x.value
        avk = ws.avk.value
        eo = ws.retrieval_eo.value
        es = ws.retrieval_ss.value
        res = ws.y.value - ws.yf.value
        #em = np.matmul(ws.dxdy.value, res)

        for rq in self.retrieval_quantities:
            rq.extract_result(x, avk, eo, es)

        return True

    def level2_diagnostics(self):
        """
        Some more advanced diagnostics tools to evaluate OEM retrievals
        Inspired from Typhon and https://github.com/maahn/pyOptimalEstimation

        https://pyoptimalestimation.readthedocs.io/en/latest/

        Work in Progress...

        """
        from typhon.retrieval.oem import error_covariance_matrix
        from scipy.linalg import inv
        ws = self.ws

        x = self.ws.x.value
        xa = self.ws.xa.value
        deltax = x-xa

        deltay = self.ws.y.value - self.ws.yf.value

        K = ws.jacobian.value
        Sx = ws.covmat_sx.value.to_dense()
        Se = ws.covmat_se.value.to_dense()
        inv_Se = np.zeros_like(Se)
        np.fill_diagonal(inv_Se, 1/np.diag(Se)) 
        self.covmat_ret =  inv(K.T @ inv_Se @ K + inv(Sx))

        KSxK = K @ Sx @ K.T
        #KSxK_inv = inv(K @ Sx @ K.T)
        #KSxK_Se = inv(KSxK + Se)

        #Syd = KSxK @ KSxK_Se @ KSxK
        
        # Testing 12.12 from Rodgers
        #chi2, chi2TestX =  _testChi2(Sx @ K.T @ KSxK_Se @ K @ Sx, deltax, significance=0.05)

        #chi2, chi2TestX =  _testChi2(Syd, deltay, significance=0.05)

        covmat_diag_ratio = np.sqrt(np.diag(self.covmat_ret) / np.diag(Sx))

        for rq in self.retrieval_quantities:

            p_grid = rq.p_grid

            avk = rq.avkm

            dof_x = np.diag(avk)
            dof = np.trace(avk)

            error_x = np.sqrt(np.diag(self.covmat_ret)[:len(p_grid)])

            error_ratio = covmat_diag_ratio[:len(p_grid)]

            plt.plot(error_ratio, np.log(p_grid))
            plt.plot(dof_x, np.log(p_grid))
            plt.gca().invert_yaxis()
            
            #self.covmat_ret = error_covariance_matrix(K,Sx,Sy )

        return True

    def get_level2_xarray(self):
        ds = xr.merge([rq.to_xarray() for rq in self.retrieval_quantities])
    
        # Spectra
        f_backend = self._sensor.f_backend if self._sensor.f_backend is not None else self.f_grid
        bad_channels = (self.noise_variance_vector > 100*np.median(self.noise_variance_vector)).astype(int)
        # measurement = self.y
        # for i, meas in enumerate(measurement):
        #     measurement[i]  = np.where(bad_channels, np.nan, meas)

        # Saving temperature profile:
        temperature = self.ws.t_field.value[:,0,0]
        original_zgrid = self.ws.z_field.value[:,0,0]
        new_zgrid = ds.o3_z.data
        interp_temperature = interpolate(new_zgrid, original_zgrid, temperature)
        
        ds['f'] = ('f', f_backend)
        ds['y'] = (('observation', 'f'), np.stack(self.y))
        ds['yf'] = (('observation', 'f'), np.stack(self.yf))
        ds['oem_diagnostics'] = ('oem_diagnostics_idx', self.oem_diagnostics)
        ds['median_noise'] = self.median_noise_level
        ds['tropospheric_opacity'] = self.tropospheric_opacity
        ds['number_of_spectra'] = self.number_of_spectra
        ds['first_sky_time'] = self.first_sky_time
        ds['last_sky_time'] = self.last_sky_time
        ds['bad_channels'] = (('observation', 'f'), np.stack(np.split(bad_channels, self.n_obs)))

        ds['temperature_profile'] = (('o3_p','o3_lat','o3_lon'), interp_temperature[:,np.newaxis,np.newaxis])

        y_baseline = self.y_baseline
        if y_baseline is not None:
            ds['y_baseline'] = (('observation', 'f'), np.stack(self.y_baseline))

        # Observations
        for field in Observation._fields:
            if field=='time':
                time_value = np.array([getattr(obs, field) for obs in self._observations])
                ds['obs_'+field] = ('observation', time_value)
            else:
                values = np.array([getattr(obs, field) for obs in self._observations], dtype=np.float)
                ds['obs_'+field] = ('observation', values)
        
        if self.n_obs < 2:
            ds = ds.expand_dims("time").assign_coords(time=("time",time_value))
            ds = ds.isel(observation=0).reset_coords('observation', drop=True)
        else: 
            # Because one can have multiple observation to a single profile.
            raise ValueError('Multiple observations geometry not covered yet')

       # ds['time'] =('observation', time_value)
       # ds = ds.swap_dims({'observation':'time'}).reset_coords('observation')
        ds.time.encoding['units'] = 'days since 2000-01-01 00:00:00'
        ds.time.encoding['calendar'] = 'proleptic_gregorian'
        # Global attributes
        ds.attrs['arts_version'] = self.arts_version
        ds.attrs['uuid'] = str(uuid.uuid4())
        ds.attrs['date_created'] = datetime.datetime.utcnow().isoformat()
        ds.attrs['nodename'] = os.uname().nodename

        return ds

    @property
    def p_grid(self):
        return self.ws.p_grid.value

    @property
    def lat_grid(self):
        return self.ws.lat_grid.value

    @property
    def lon_grid(self):
        return self.ws.lon_grid.value

    @property
    def f_grid(self):
        return self.ws.f_grid.value

    @property
    def n_p(self):
        return len(self.ws.p_grid.value)

    @property
    def n_lat(self):
        if self.atmosphere_dim == 1:
            return 1
        return len(self.ws.lat_grid.value)

    @property
    def n_lon(self):
        if self.atmosphere_dim == 1:
            return 1
        return len(self.ws.lat_grid.value)

    @property
    def n_y(self):
        return len(self.ws.y.value)

    @property
    def y(self):
        y = np.copy(self.ws.y.value)
        return np.split(y, self.n_obs)

    @property
    def yf(self):
        yf = np.copy(self.ws.yf.value)
        return np.split(yf, self.n_obs)

    @property
    def median_noise_level(self):
        return np.sqrt(np.median(self.noise_variance_vector))

    @property
    def y_baseline(self):
        if self.ws.y_baseline.initialized:
            bl = np.copy(self.ws.y_baseline.value)
            if bl.size == 1 and bl[0] == 0:
                return None
            return np.split(bl, self.n_obs)
        else:
            return None

    @property
    def n_obs(self):
        if self.ws.sensor_time.initialized:
            return len(self.ws.sensor_time.value)
        else:
            return len(self.ws.sensor_los.value)

    @property
    def abs_species_maintags(self):
        abs_species = self.ws.abs_species.value
        maintags = [str(st[0]).split('-')[0] for st in abs_species]
        return maintags

    @property
    def atmosphere_dim(self):
        return self.ws.atmosphere_dim.value

    @property
    def oem_converged(self):
        return self.ws.oem_diagnostics.value[0] == 0

    @property
    def oem_diagnostics(self):
        return self.ws.oem_diagnostics.value

    @property
    def oem_errors(self):
        return self.ws.oem_errors.value

    @property
    def arts_version(self):
        cmd = os.path.join(os.environ['ARTS_BUILD_PATH'], 'src/arts')
        output = os.popen(cmd + ' --version').read().splitlines()
        return output[0]

def _estimateChi2(S, z, atol=1e-5):
    '''
    From : https://github.com/maahn/pyOptimalEstimation/blob/master/pyOptimalEstimation/pyOEcore.py
    Estimate Chi^2 to estimate whether z is from distribution with 
    covariance S
    Parameters
    ----------
    S : {array}
        Covariance matrix
    z : {array}
        Vector to test
        atol : float (default 1e-5)
            The absolute tolerance for comparing eigen values to zero. We 
            found that values should be than the numpy.isclose defualt value 
            of 1e-8.
    Returns
    -------
    float
        Estimated chi2 value
    '''

    eigVals, eigVecsL = scipy.linalg.eig(S, left=True, right=False)
    z_prime = eigVecsL.T.dot(z)

    # Handle singular matrices. See Rodgers ch 12.2
    notNull = np.abs(eigVals) > atol
    dofs = np.sum(notNull)
    if dofs != len(notNull):
        print('Warning. Singular Matrix with rank %i instead of %i. '
              '(This is typically save to ignore)       ' %
              (dofs, len(notNull)))

    # Rodgers eq. 12.1
    chi2s = z_prime[notNull]**2/eigVals[notNull]
    return chi2s, dofs


def _testChi2(S, z, significance, atol=1e-5):
    '''
    From : https://github.com/maahn/pyOptimalEstimation/blob/master/pyOptimalEstimation/pyOEcore.py
    Test whether z is from distribution with covariance S with significance
    Parameters
    ----------
    S : {array}
        Covariance matrix
    z : {array}
        Vector to test
    significance : {float}
        Significane level
        atol : float (default 1e-5)
            The absolute tolerance for comparing eigen values to zero. We 
            found that values should be than the numpy.isclose defualt value 
            of 1e-8.
    Returns
    -------
    float
        Estimated chi2 value
    float
        Theoretical chi2 value for significance
    bool
        True if Chi^2 test passed
    '''
    chi2s_obs, dof = _estimateChi2(S, z, atol=atol)
    chi2_obs = np.real_if_close(np.sum(chi2s_obs))
    chi2_theo = scipy.stats.chi2.isf(significance, dof)
    # chi2_theo1 = scipy.stats.chi2.isf(significance, 1)

    # print(chi2_obs<= chi2_theo, np.all(chi2s_obs<= chi2_theo1))

    return chi2_obs, chi2_theo