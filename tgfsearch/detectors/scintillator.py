"""A class for keeping track of a particular scintillator's data."""
import pandas as pd
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
    frame : pd.DataFrame
        A pandas dataframe containing all the scintillator's data.
    calibration : list
        A list containing calibration info. The energy bins corresponding to compton shoulders/photo peaks
    filelist : list
        A list of scintillator files for a particular day.
    filetime_extrema : list
        A list of lists. Each sublist contains a pair of numbers corresponding to the first and last second in each file
    passtime : dict
        A dictionary containing information needed to import the subsequent files properly (if applicable).

    """

    def __init__(self, name, eRC):
        self.name = name
        self.eRC = eRC
        self.frame = pd.DataFrame()
        self.calibration = []
        self.filelist = []
        self.filetime_extrema = []
        self.passtime = {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0, 'ppsunix': -1.0, 'lastwc': 0,
                         'ppswc': 0, 'hz': 8e7, 'started': 0}

    def __str__(self):
        return f'Scintillator({self.name}, {self.eRC})'

    def __repr__(self):
        return self.__str__()

    def __bool__(self):
        # Using energy as an arbitrary check here. Time or wc would've worked too.
        return True if 'energy' in self.frame and len(self.frame['energy']) > 0 else False

    def get_attribute(self, attribute):
        if attribute in self.frame:
            return self.frame[attribute].to_numpy()
        elif hasattr(self, attribute):
            return getattr(self, attribute)
        else:
            raise ValueError(f"'{attribute}' is either not a valid attribute or data hasn't been imported.")

    def set_attribute(self, attribute, info):
        if attribute in self.frame:
            info_type = type(info)
            if info_type == np.ndarray or info_type == list:
                if len(info) == len(self.frame[attribute]):
                    self.frame[attribute] = info
                else:
                    raise ValueError(f'length of info ({len(info)}) does not match the number '
                                     f'of scintillator frame indices ({len(self.frame[attribute])}).')

            else:
                raise TypeError(f"'{attribute}' must be a numpy array or list, not '{info_type.__name__}'.")

        elif hasattr(self, attribute):
            attribute_type = type(getattr(self, attribute))
            info_type = type(info)
            if info_type == attribute_type:
                setattr(self, attribute, info)
            else:
                raise TypeError(f"'{attribute}' must be of type '{attribute_type.__name__}', "
                                f"not '{info_type.__name__}'.")

        else:
            raise ValueError(f"'{attribute}' is either not a valid attribute or data hasn't been imported.")
