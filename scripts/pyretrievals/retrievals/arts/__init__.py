"""
This module implements some classes and helpers to deal with ARTS in a more pythonic
way.
"""

from dotenv import load_dotenv


load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')


# Interface stuff
from .interface import (
    Observation,
    ArtsController,
    OemException,
)


# Atmosphere
from .atmosphere import (
    Atmosphere,
)


# Sensors
from .sensors import (
    SensorOff,
    SensorGaussian,
    SensorFFT,
    SensorFFT_Sideband
)


# Retrievals
from .retrieval import (
    RetrievalQuantity,
    GriddedRetrievalQuantity,
    AbsSpecies,
    Wind,
    FreqShift,
    Polyfit,
    Sinefit
)
