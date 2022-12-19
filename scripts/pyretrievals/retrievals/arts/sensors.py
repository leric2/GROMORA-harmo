from abc import ABC, abstractmethod, abstractstaticmethod
import numpy as np

#from typhon.arts.workspace import arts_agenda
#from typhon.arts.griddedfield import GriddedField1
from pyarts.workspace import arts_agenda
from pyarts.griddedfield import GriddedField1
import matplotlib.pyplot as plt

class AbstractSensor(ABC):
    """
    Abstract Sensor requires the implementation of an `sensor_response_agenda` property.
    See the specific derived classes for how to use a sensor.
    """

    def apply(self, ws):
        """Copy and execute sensor response agenda."""
        ws.Copy(ws.sensor_response_agenda, self.sensor_response_agenda)
        ws.AgendaExecute(ws.sensor_response_agenda)

    @property
    @abstractmethod
    def sensor_response_agenda(self):
        """
        The sensor response agenda.

        :type: typhon.arts.workspace.agendas.Agenda

        .. seealso:: Decorator :py:func:`typhon.arts.workspace.workspace.arts_agenda`.
        """
        pass

    @property
    @abstractmethod
    def f_backend(self):
        pass


class SensorOff(AbstractSensor):
    """Sensor that does nothing."""

    def __init__(self):
        pass

    def apply(self, ws):
        """Copy and execute sensor response agenda."""
        ws.AntennaOff()
        ws.sensorOff()
        pass

    @property
    def sensor_response_agenda(self):
        return None

    @property
    def f_backend(self):
        return None


class SensorFFT(AbstractSensor):
    """
    Sensor with channel response for an FFT Spectrometer with :math:`\mathrm{sinc}^2` response.
    """

    def __init__(self, f_backend, resolution, num_channels=10):
        """
        :param f_backend: The backend frequency vector.
        :param resolution: The frequency resolution of the FFTS in Hz
        :param num_channels: Number of channels with nonzero response, default: 10
        """
        self._f_backend = f_backend
        self.resolution = resolution
        self.num_channels = num_channels

        # Compute the backend channel response
        grid = np.linspace(-self.num_channels / 2,
                           self.num_channels / 2,
                           20 * self.num_channels)
        response = np.sinc(grid) ** 2
        self.bcr = GriddedField1(name='Backend channel response function for FFTS',
                                 gridnames=['Frequency'], dataname='Data',
                                 grids=[self.resolution * grid],
                                 data=response)

    def apply(self, ws):
        # Modify workspace
        ws.FlagOn(ws.sensor_norm)
        ws.f_backend = self.f_backend
        ws.backend_channel_response = [self.bcr, ]

        super().apply(ws)

    @property
    def sensor_response_agenda(self):
        @arts_agenda
        def sensor_response_agenda(ws):
            ws.AntennaOff()
            ws.sensor_responseInit()
            ws.sensor_responseBackend()

        return sensor_response_agenda

    @property
    def f_backend(self):
        return self._f_backend

class SensorFFT_Sideband(AbstractSensor):
    """
    Sensor with channel response for an FFT Spectrometer with :math:`\mathrm{sinc}^2` response.

    In addition, this sensor include a sideband response. 

    """

    def __init__(self, f_backend, resolution, num_channels, lo_freq, sideband_mode, intermediate_freq, sideband_response):
        """
        :param f_backend: The backend frequency vector.
        :param resolution: The frequency resolution of the FFTS in Hz
        :param num_channels: Number of channels with nonzero response, default: 10
        :param lo_freq: The local oscillator frequency in Hz
        :param sideband_mode: the type of sideband, upper or lower
        :param intermediate_freq: the intermediate frequency of the sideband
        :param sideband_response: the sideband response corresponding to the intermediate frequency
        """
        self._f_backend = f_backend
        self.resolution = resolution
        self.num_channels = num_channels
        self.lo_freq = lo_freq
        self.sideband_mode = sideband_mode

        # Compute the backend channel response
        grid = np.linspace(-self.num_channels / 2,
                           self.num_channels / 2,
                           20 * self.num_channels)
        response = np.sinc(grid) ** 2
        self.bcr = GriddedField1(name='Backend channel response function for FFTS',
                                 gridnames=['Frequency'], dataname='Data',
                                 grids=[self.resolution * grid],
                                 data=response)

        self.sideband_response = GriddedField1(name='Sideband response function for mixer',
                         gridnames=['Frequency'], dataname='Data',
                         grids=[intermediate_freq],
                          data=sideband_response)


    def apply(self, ws):
        # Modify workspace
        ws.FlagOn(ws.sensor_norm)
        ws.f_backend = self.f_backend
        ws.backend_channel_response = [self.bcr, ]
        ws.lo = self.lo_freq
        ws.sideband_mode = self.sideband_mode
        ws.sideband_response = self.sideband_response

        super().apply(ws)

    @property
    def sensor_response_agenda(self):
        @arts_agenda
        def sensor_response_agenda(ws):
            ws.AntennaOff()
            ws.sensor_responseInit()
            ws.sensor_responseMixer(lo=self.lo_freq)
            ws.sensor_responseIF2RF(sideband_mode=self.sideband_mode)
            ws.sensor_responseBackend()

        return sensor_response_agenda

    @property
    def f_backend(self):
        return self._f_backend


class SensorFFT_Sideband_Antenna(AbstractSensor):
    """
    Sensor with channel response for an FFT Spectrometer with :math:`\mathrm{sinc}^2` response.

    In addition, this sensor include a sideband response. 

    """

    def __init__(self, f_backend, resolution, num_channels, lo_freq, sideband_mode, intermediate_freq, sideband_response, fwhm):
        """
        :param f_backend: The backend frequency vector.
        :param resolution: The frequency resolution of the FFTS in Hz
        :param num_channels: Number of channels with nonzero response, default: 10
        :param lo_freq: The local oscillator frequency in Hz
        :param sideband_mode: the type of sideband, upper or lower
        :param intermediate_freq: the intermediate frequency of the sideband
        :param sideband_response: the sideband response corresponding to the intermediate frequency
        :param fwhm: the full width half maximum of the antenna to consider (in degree)
        """
        self._f_backend = f_backend
        self.resolution = resolution
        self.num_channels = num_channels
        self.lo_freq = lo_freq
        self.sideband_mode = sideband_mode
        self.fmhw = fwhm

        # Compute the backend channel response
        grid = np.linspace(-self.num_channels / 2,
                           self.num_channels / 2,
                           20 * self.num_channels)
        response = np.sinc(grid) ** 2
        self.bcr = GriddedField1(name='Backend channel response function for FFTS',
                                 gridnames=['Frequency'], dataname='Data',
                                 grids=[self.resolution * grid],
                                 data=response)

        self.sideband_response = GriddedField1(name='Sideband response function for mixer',
                         gridnames=['Frequency'], dataname='Data',
                         grids=[intermediate_freq],
                          data=sideband_response)


    def apply(self, ws):
        # Modify workspace
        ws.FlagOn(ws.sensor_norm)
        ws.f_backend = self.f_backend
        ws.backend_channel_response = [self.bcr, ]
        ws.lo = self.lo_freq
        ws.sideband_mode = self.sideband_mode
        ws.sideband_response = self.sideband_response
        #ws.AntennaConstantGaussian1D(n_za_grid=2, fwhm=self.fmhw, xwidth_si=3, dx_si=0.2)
        #
        # ws.ReadXML(ws.antenna_response, "/home/esauvageat/Documents/GROMORA/Analysis/GROMORA-harmo/files/antenna_response.xml")
        # ws.mblock_dlos_gridUniformCircular(spacing=0.2, width=3)
        # ws.antenna_dlos = [0]
        # ws.antenna_dim =  1
        # ws.mblock_dlos_grid = np.array([-2,2])
        
        super().apply(ws)

    @property
    def sensor_response_agenda(self):
        @arts_agenda
        def sensor_response_agenda(ws):
            ws.AntennaConstantGaussian1D(n_za_grid=2, fwhm=self.fmhw, xwidth_si=3, dx_si=0.5)
            #ws.antenna_responseGaussian(fwhm=self.fmhw, xwidth_si=3, dx_si=0.5)
            #ws.mblock_dlos_grid = np.array([-2,2])
            #ws.AntennaConstantGaussian1D(n_za_grid=1, fwhm=1.5)
            #print(ws.mblock_dlos_grid)
            #ws.AntennaMultiBeamsToPencilBeams()
            # ws.ReadXML(ws.antenna_response, "/home/esauvageat/Documents/GROMORA/Analysis/GROMORA-harmo/files/antenna_response.xml")
            #ws.mblock_dlos_gridUniformRectangular(spacing=0.2, za_width=0.1, aa_width=0)
            #ws.antenna_dlos = [0]
            ws.sensor_responseInit()
            ws.sensor_responseAntenna()
            ws.sensor_responseMixer(lo=self.lo_freq)
            ws.sensor_responseIF2RF(sideband_mode=self.sideband_mode)
            ws.sensor_responseBackend()

        return sensor_response_agenda

    @property
    def f_backend(self):
        return self._f_backend


class SensorGaussian(AbstractSensor):
    """Sensor with Gaussian Channel response."""

    def __init__(self, f_backend, fwhm):
        """
        :param f_backend: Backend frequencies
        :param fwhm: Full width at half maximum (resolution)
        """
        self._f_backend = f_backend
        self.fwhm = fwhm

    def apply(self, ws):
        ws.FlagOn(ws.sensor_norm)
        ws.f_backend = self.f_backend
        ws.backend_channel_responseGaussian(fwhm=self.fwhm) #xwidth_si=0.5
        super().apply(ws)

    @property
    def sensor_response_agenda(self):
        @arts_agenda
        def sensor_response_agenda(ws):
            ws.AntennaOff()
            ws.sensor_responseInit()
            ws.sensor_responseBackend()

        return sensor_response_agenda

    @property
    def f_backend(self):
        return self._f_backend

class SensorRectSB_Antenna(AbstractSensor):
    """Sensor with rectangular Channel response, now including gaussian Antenna patter. Used for FB retrievals"""

    def __init__(self, f_backend, channel_width, lo_freq, sideband_mode, intermediate_freq, sideband_response, fwhm):
        """
        :param f_backend: Backend frequencies
        """
        # Compute the backend channel response
        self._f_backend = f_backend
        self._channel_width = channel_width
        self.lo_freq = lo_freq
        self.sideband_mode = sideband_mode
        self.bcr = []
        self.fmhw = fwhm
        for i, df in enumerate(channel_width):
            self.bcr.append(
                GriddedField1(
                    name='Backend channel response function for FB',
                    gridnames=['Frequency'], 
                    dataname='Data',
                    grids=[np.array((-df/2 , df/2))],
                    data=np.array((1,1))))
        
        self.sideband_response= GriddedField1(
                    name='Sideband response function for mixer',
                    gridnames=['Frequency'], 
                    dataname='Data',
                    grids=[intermediate_freq],
                    data=sideband_response)

    def apply(self, ws):
        ws.FlagOn(ws.sensor_norm)
        ws.f_backend = self.f_backend
        ws.backend_channel_response = self.bcr
        ws.lo = self.lo_freq
        ws.sideband_mode = self.sideband_mode
        ws.sideband_response = self.sideband_response
        
        super().apply(ws)

    @property
    def sensor_response_agenda(self):
        @arts_agenda
        def sensor_response_agenda(ws):
            ws.AntennaConstantGaussian1D(n_za_grid=2, fwhm=self.fmhw, xwidth_si=3, dx_si=0.5)
            ws.sensor_responseInit()
            ws.sensor_responseAntenna()
            ws.sensor_responseMixer(lo=self.lo_freq)
            ws.sensor_responseIF2RF(sideband_mode=self.sideband_mode)
            ws.sensor_responseBackend()

        return sensor_response_agenda


    @property
    def f_backend(self):
        return self._f_backend

class SensorRectSB(AbstractSensor):
    """Sensor with rectangular Channel response. Used for FB retrievals"""

    def __init__(self, f_backend, channel_width, lo_freq, sideband_mode, intermediate_freq, sideband_response):
        """
        :param f_backend: Backend frequencies
        """
        # Compute the backend channel response
        self._f_backend = f_backend
        self._channel_width = channel_width
        self.lo_freq = lo_freq
        self.sideband_mode = sideband_mode
        self.bcr = []
        for i, df in enumerate(channel_width):
            self.bcr.append(
                GriddedField1(
                    name='Backend channel response function for FB',
                    gridnames=['Frequency'], 
                    dataname='Data',
                    grids=[np.array((-df/2 , df/2))],
                    data=np.array((1,1))))
        
        self.sideband_response= GriddedField1(
                    name='Sideband response function for mixer',
                    gridnames=['Frequency'], 
                    dataname='Data',
                    grids=[intermediate_freq],
                    data=sideband_response)

    def apply(self, ws):
        ws.FlagOn(ws.sensor_norm)
        ws.f_backend = self.f_backend
        ws.backend_channel_response = self.bcr
        ws.lo = self.lo_freq
        ws.sideband_mode = self.sideband_mode
        ws.sideband_response = self.sideband_response
        
        super().apply(ws)

    @property
    def sensor_response_agenda(self):
        @arts_agenda
        def sensor_response_agenda(ws):
            ws.AntennaOff()
            ws.sensor_responseInit()
            ws.sensor_responseMixer(lo=self.lo_freq)
            ws.sensor_responseIF2RF(sideband_mode=self.sideband_mode)
            ws.sensor_responseBackend()

        return sensor_response_agenda


    @property
    def f_backend(self):
        return self._f_backend

class SensorRect(AbstractSensor):
    """Sensor with rectangular Channel response. Used for FB retrievals"""

    def __init__(self, f_backend, channel_width):
        """
        :param f_backend: Backend frequencies
        """
        # Compute the backend channel response
        self._f_backend = f_backend
        self._channel_width = channel_width
        self.bcr = []
        for i, df in enumerate(channel_width):
             self.bcr.append(GriddedField1(name='Backend channel response function for FB',
                                 gridnames=['Frequency'], dataname='Data',
                                 grids=[np.array((-df/2 , df/2))],
                                 data=np.array((1,1))))
        
    def apply(self, ws):
        ws.FlagOn(ws.sensor_norm)
        ws.f_backend = self.f_backend
        ws.backend_channel_response = self.bcr
        
        super().apply(ws)

    @property
    def sensor_response_agenda(self):
        @arts_agenda
        def sensor_response_agenda(ws):
            ws.AntennaOff()
            ws.sensor_responseInit()
            ws.sensor_responseBackend()

        return sensor_response_agenda

    @property
    def f_backend(self):
        return self._f_backend