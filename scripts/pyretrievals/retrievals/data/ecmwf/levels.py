import numpy as np
import xarray as xr


def pressure_levels(lnsp):
    """
    If vector `lnsp` has dims (time,) then return an array with dims (level, time) containing pressure.
    """
    a = hybrid_level_a[:, np.newaxis]
    b = hybrid_level_b[:, np.newaxis]
    sp = np.exp(lnsp[np.newaxis,:])  # surface pressure
    return a + b * sp

def pressure_levels_91(lnsp):
    """
    If vector `lnsp` has dims (time,) then return an array with dims (level, time) containing pressure.
    """
    a = hybrid_level_a_91[:, np.newaxis]
    b = hybrid_level_b_91[:, np.newaxis]
    sp = np.exp(lnsp[np.newaxis,:])  # surface pressure
    return a + b * sp


def level_to_pressure(level, lnsp=0):
    """
    Return the pressure for a given level and logarithm of surface pressure `lnsp`.
    If hybrid level `b` parameter, `lnsp` is ignored.
    """
    a = hybrid_level_a[level-1]
    b = hybrid_level_b[level-1]
    if b != 0 and lnsp == 0:
        raise ValueError('Hybrid level b parameter is not zero.')
    return a + b * np.exp(lnsp)

def level_to_pressure_91(level, lnsp=0):
    """
    Return the pressure for a given level and logarithm of surface pressure `lnsp`.
    If hybrid level `b` parameter, `lnsp` is ignored.
    """
    a = hybrid_level_a_91[level-1]
    b = hybrid_level_b_91[level-1]
    if b != 0 and lnsp == 0:
        raise ValueError('Hybrid level b parameter is not zero.')
    return a + b * np.exp(lnsp)


""" ECMWF hybrid level definition <https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels> """
hybrid_level_n = np.arange(1, 138)
hybrid_level_n_91 = np.arange(1, 92)


#: ECMWF hybrid level parameter `a`
hybrid_level_a = np.array(
    [2.000365, 3.102241, 4.666084, 6.827977, 9.746966, 13.605424, 18.608931, 24.985718, 32.98571, 42.879242,
     54.955463, 69.520576, 86.895882, 107.415741, 131.425507, 159.279404, 191.338562, 227.968948, 269.539581,
     316.420746, 368.982361, 427.592499, 492.616028, 564.413452, 643.339905, 729.744141, 823.967834, 926.34491,
     1037.201172, 1156.853638, 1285.610352, 1423.770142, 1571.622925, 1729.448975, 1897.519287, 2076.095947,
     2265.431641, 2465.770508, 2677.348145, 2900.391357, 3135.119385, 3381.743652, 3640.468262, 3911.490479,
     4194.930664, 4490.817383, 4799.149414, 5119.89502, 5452.990723, 5798.344727, 6156.074219, 6526.946777,
     6911.870605, 7311.869141, 7727.412109, 8159.354004, 8608.525391, 9076.400391, 9562.682617, 10065.978516,
     10584.631836, 11116.662109, 11660.067383, 12211.547852, 12766.873047, 13324.668945, 13881.331055,
     14432.139648, 14975.615234, 15508.256836, 16026.115234, 16527.322266, 17008.789063, 17467.613281,
     17901.621094, 18308.433594, 18685.71875, 19031.289063, 19343.511719, 19620.042969, 19859.390625,
     20059.931641, 20219.664063, 20337.863281, 20412.308594, 20442.078125, 20425.71875, 20361.816406,
     20249.511719, 20087.085938, 19874.025391, 19608.572266, 19290.226563, 18917.460938, 18489.707031,
     18006.925781, 17471.839844, 16888.6875, 16262.046875, 15596.695313, 14898.453125, 14173.324219,
     13427.769531, 12668.257813, 11901.339844, 11133.304688, 10370.175781, 9617.515625, 8880.453125, 8163.375,
     7470.34375, 6804.421875, 6168.53125, 5564.382813, 4993.796875, 4457.375, 3955.960938, 3489.234375,
     3057.265625, 2659.140625, 2294.242188, 1961.5, 1659.476563, 1387.546875, 1143.25, 926.507813, 734.992188,
     568.0625, 424.414063, 302.476563, 202.484375, 122.101563, 62.78125, 22.835938, 3.757813, 0, 0])

# ECMWF before 25.06.2013
hybrid_level_a_91 = np.array(
    [2.000040, 3.980832, 7.387186, 12.908319, 21.413612,33.952858,51.746601,76.167656,108.715561,150.986023,204.637451,
    271.356506,352.824493,450.685791,566.519226,701.813354,857.945801,1036.166504,1237.585449,1463.163940,1713.709595,1989.874390,
    2292.155518,2620.898438,2976.302246,3358.425781,3767.196045,4202.416504,4663.776367,5150.859863,5663.156250,6199.839355,6759.727051,
    7341.469727,7942.926270,8564.624023,9208.305664,9873.560547,10558.881836,11262.484375,11982.662109,12713.897461,13453.225586,
    14192.009766,14922.685547,15638.053711,16329.560547,16990.623047,17613.281250,18191.029297, 18716.96875,19184.544922,19587.513672,19919.796875,
    20175.394531,20348.916016,20434.158203,20426.21875,20319.011719,20107.03125,19785.357422,19348.775391,18798.822266,18141.296875,17385.595703,
    16544.585938,15633.566406,14665.645508,13653.219727,12608.383789,11543.166992,10471.310547,9405.222656,8356.25293,7335.164551,6353.920898,
    5422.802734,4550.21582,3743.464355,3010.146973,2356.202637,1784.854614,1297.656128,895.193542,576.314148,336.772369,162.043427,54.208336,
    6.575628,0.00316,0]
)

#: ECMWF hybrid level parameter `b`
hybrid_level_b = np.array(
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7E-6, 2.4E-5, 5.9E-5, 0.000112, 0.000199, 0.00034, 0.000562,
     0.00089, 0.001353, 0.001992, 0.002857, 0.003971, 0.005378, 0.007133, 0.009261, 0.011806, 0.014816, 0.018318,
     0.022355, 0.026964, 0.032176, 0.038026, 0.044548, 0.051773, 0.059728, 0.068448, 0.077958, 0.088286, 0.099462,
     0.111505, 0.124448, 0.138313, 0.153125, 0.16891, 0.185689, 0.203491, 0.222333, 0.242244, 0.263242, 0.285354,
     0.308598, 0.332939, 0.358254, 0.384363, 0.411125, 0.438391, 0.466003, 0.4938, 0.521619, 0.549301, 0.576692,
     0.603648, 0.630036, 0.655736, 0.680643, 0.704669, 0.727739, 0.749797, 0.770798, 0.790717, 0.809536, 0.827256,
     0.843881, 0.859432, 0.873929, 0.887408, 0.8999, 0.911448, 0.922096, 0.931881, 0.94086, 0.949064, 0.95655, 0.963352,
     0.969513, 0.975078, 0.980072, 0.984542, 0.9885, 0.991984, 0.995003, 0.99763, 1])


# ECMWF before 25.06.2013
hybrid_level_b_91 = np.array(
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000014,0.000055,0.000131,0.000279,0.000548,0.001,
    0.001701,0.002765,0.004267,0.006322,0.009035,0.012508,0.01686,0.022189,0.02861,0.036227,0.045146,0.055474,0.067316,0.080777,
    0.095964,0.112979,0.131935,0.152934,0.176091,0.20152,0.229315,0.259554,0.291993,0.326329,0.362203,0.399205,0.436906,0.475016,
    0.51328,0.551458,0.589317,0.626559,0.662934,0.698224,0.732224,0.764679,0.795385,0.824185,0.85095,0.875518,0.897767,0.917651,
    0.935157,0.950274,0.963007,0.973466,0.982238,0.989153,0.994204,0.99763,1]
)


#: ECMWF hybrid level parameters
hybrid_level = xr.Dataset({'a': ('level', hybrid_level_a), 'b': ('level', hybrid_level_b)}, coords={'level': hybrid_level_n})
hybrid_level_91 = xr.Dataset({'a': ('level', hybrid_level_a_91), 'b': ('level', hybrid_level_b_91)}, coords={'level': hybrid_level_n_91})