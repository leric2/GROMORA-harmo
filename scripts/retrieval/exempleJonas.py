#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 10:17:03 2020

@author: eric
"""
from abc import ABC


class DataLevel0(ABC):
    def __init__(self, filename=None):
        self._filename = filename

        # some dummy data
        self._data = dict()
        self._data["Temperature.0"] = 10
        self._data["Temperature.1"] = 20
        self._data["Temperature.2"] = 30

    def get_temerature_reading(self, time, channel):
        """ Get the temperature form a specific temperature channel.
        """
        return self._data[f"Temperature.{channel}"]

    def get_hot_load_temperature(self, time):
        """ Get hot load temperature for a specific time.
        """
        raise NotImplementedError("abstract base class")

class GROMOS_Level0(DataLevel0):
    def get_hot_load_temperature(self, time):
        """ On GROMOS, hot load temperature is on channel 0.
        """
        return self.get_temerature_reading(time, 0)

class SOMORA_Level0(DataLevel0):
    def get_hot_load_temperature(self, time):
        """ On SOMORA, hot load temperature is on channel 2.
        """
        return self.get_temerature_reading(time, 2)


def evaluate(instrument):
    # This function does something fancy
    hltemp = instrument.get_hot_load_temperature(1234556)
    return hltemp


if __name__ == "__main__":
    l0_gromos = GROMOS_Level0("filename.txt")
    l0_somora = SOMORA_Level0("filename.txt")

    print("GROMOS hot load:", evaluate(l0_gromos))
    print("SOMORA hot load:", evaluate(l0_somora))

