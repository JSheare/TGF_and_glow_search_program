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
    lm_frame : pd.DataFrame
        A pandas dataframe containing all the scintillator's data.
    calibration_bins : list
        A list containing the energy bins corresponding to Compton edges/photo peaks used for calibration.
    calibration_energies : list
        A list containing the energies of Compton edges/photo peaks used for calibration.
    lm_filelist : list
        A list of list mode files for the day.
    lm_filetime_extrema : list
        A list of lists. Each sublist contains a pair of numbers corresponding to
        the first and last second in each list mode file.
    trace_filelist : list
        A list of trace files for the day.
    traces : dict
        A dictionary containing trace data for each of the day's trace files.
    passtime : dict
        A dictionary containing information needed to import the subsequent files properly (if applicable).

    """

    def __init__(self, name, eRC):
        self.name = name
        self.eRC = eRC
        self.lm_frame = pd.DataFrame()
        self.calibration_energies = []
        self.calibration_bins = []
        self.lm_filelist = []
        self.lm_filetime_extrema = []
        self.trace_filelist = []
        self.traces = {}
        self.passtime = {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0, 'ppsunix': -1.0, 'lastwc': 0,
                         'ppswc': 0, 'hz': 8e7, 'started': 0}

    def __str__(self):
        return f'Scintillator({self.name}, {self.eRC})'

    def __repr__(self):
        return self.__str__()

    def __bool__(self):
        # Using energy as an arbitrary check here. Time or wc would've worked too.
        return True if 'energy' in self.lm_frame and len(self.lm_frame['energy']) > 0 else False

    def get_attribute(self, attribute):
        if attribute in self.lm_frame:
            return self.lm_frame[attribute].to_numpy()
        elif hasattr(self, attribute):
            return getattr(self, attribute)
        else:
            raise ValueError(f"'{attribute}' is either not a valid attribute or data hasn't been imported.")

    def set_attribute(self, attribute, info):
        if attribute in self.lm_frame:
            info_type = type(info)
            if info_type == np.ndarray or info_type == list:
                if len(info) == len(self.lm_frame[attribute]):
                    self.lm_frame[attribute] = info
                else:
                    raise ValueError(f'length of info ({len(info)}) does not match the number '
                                     f'of scintillator frame indices ({len(self.lm_frame[attribute])}).')

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

    def get_trace_list(self):
        return list(self.traces.keys())

    def get_trace(self, time_id):
        if time_id in self.traces:
            return self.traces[time_id]
        else:
            raise ValueError(f"No trace with name '{time_id}' for scintillator '{self.name}'.")
