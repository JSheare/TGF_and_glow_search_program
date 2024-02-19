"""A class for keeping track of a particular scintillator's data."""
import numpy as np


class Scintillator:
    """Used to store data for a single scintillator.

    Parameters
    ----------
    name : str
        The scintillator's name (usually an abbreviation).
    eRC : str
        The scintillator's serial number.

    Attributes
    ----------
    calibration : list
        A list containing calibration info. The energy bins corresponding to compton shoulders/photo peaks
    filelist : list
        A list of scintillator files for a particular day.
    filetime_extrema : list
        A list of lists. Each sublist contains a pair of numbers corresponding to the first and last second in each file
    time : np.array
        Contains second-of-day time data.
    energy : np.array
        Contains energy data.
    wc : np.array
        Contains wallclock time data.
    passtime : dict
        A dictionary containing information needed to import the subsequent files properly (if applicable).

    """

    def __init__(self, name, eRC):
        self.name = name
        self.eRC = eRC
        self.calibration = []
        self.filelist = []
        self.filetime_extrema = []
        self.time = np.array([])
        self.energy = np.array([])
        self.wc = np.array([])
        self.passtime = {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0, 'ppsunix': -1.0, 'lastwc': 0,
                         'ppswc': 0, 'hz': 8e7, 'started': 0}

    def __str__(self):
        return f'Scintillator({self.name}, {self.eRC})'

    def __repr__(self):
        return self.__str__()

    def __bool__(self):
        # Using energy as an arbitrary check here. Time or wc would've worked too.
        return True if len(self.energy) > 0 else False

    def get_attribute(self, attribute):
        if hasattr(self, attribute):
            return getattr(self, attribute)
        else:
            raise ValueError(f"'{attribute}' is not a valid attribute.")

    def set_attribute(self, attribute, info):
        if hasattr(self, attribute):
            setattr(self, attribute, info)
        else:
            raise ValueError(f"'{attribute}' is not a valid attribute.")
