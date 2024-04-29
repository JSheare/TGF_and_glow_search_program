"""A class for keeping track of a particular scintillator's data."""
import copy as copy
import pandas as pd
import numpy as np

import tgfsearch.parameters as params
import tgfsearch.tools as tl


class Scintillator:
    """A class used to store data for a single scintillator.

    Parameters
    ----------
    name : str
        The scintillator's name (abbreviated).
    eRC : str
        The scintillator's serial number.

    Attributes
    ----------
    lm_frame : pandas.core.frame.Dataframe
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
    lm_file_indices : dict
        A dictionary of lists. Each list contains the indices needed to slice data for a particular file
        out of lm_frame.
    trace_filelist : list
        A list of trace files for the day.
    traces : dict
        A dictionary containing trace data for each of the day's trace files.
    passtime : dict
        A dictionary containing information needed to import the subsequent list mode files properly (if applicable).

    """

    def __init__(self, name, eRC):
        self.name = name
        self.eRC = eRC
        self.lm_frame = pd.DataFrame()
        self.calibration_energies = []
        self.calibration_bins = []
        self.lm_filelist = []
        self.lm_filetime_extrema = []
        self.lm_file_indices = {}
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
        """Returns the requested attribute."""
        if attribute in self.lm_frame:
            return self.lm_frame[attribute].to_numpy()
        elif hasattr(self, attribute):
            return copy.deepcopy(getattr(self, attribute))

        else:
            raise ValueError(f"'{attribute}' is either not a valid attribute or data hasn't been imported.")

    def set_attribute(self, attribute, info):
        """Updates the requested attribute."""
        if attribute in self.lm_frame:
            info_type = type(info)
            if info_type == np.ndarray or info_type == list:
                if len(info) == len(self.lm_frame[attribute]):
                    self.lm_frame[attribute] = copy.deepcopy(info)
                else:
                    raise ValueError(f'length of info ({len(info)}) does not match the number '
                                     f'of scintillator frame indices ({len(self.lm_frame[attribute])}).')

            else:
                raise TypeError(f"'{attribute}' must be a numpy array or list, not '{info_type.__name__}'.")

        elif hasattr(self, attribute):
            attribute_type = type(getattr(self, attribute))
            info_type = type(info)
            if info_type == attribute_type:
                if attribute == 'lm_filelist':
                    info = tl.filter_files(info)  # To ensure that find_lm_file_index works properly

                setattr(self, attribute, copy.deepcopy(info))
            else:
                raise TypeError(f"'{attribute}' must be of type '{attribute_type.__name__}', "
                                f"not '{info_type.__name__}'.")

        else:
            raise ValueError(f"'{attribute}' is either not a valid attribute or data hasn't been imported.")

    def find_lm_file_index(self, count_time):
        """Returns the index of the list mode file that the given count occurred in."""
        # Checking to see that the count is inside the day or in the ~500 seconds of the next day included sometimes
        if count_time < 0 or count_time > params.SEC_PER_DAY + 500:
            return -1

        # Binary search of list mode file ranges
        low = 0
        high = len(self.lm_filetime_extrema) - 1
        while low <= high:
            mid = low + (high - low) // 2
            if self.lm_filetime_extrema[mid][0] <= count_time <= self.lm_filetime_extrema[mid][1]:
                return mid
            elif self.lm_filetime_extrema[mid][0] > count_time:
                high = mid - 1
            else:
                low = mid + 1

        return -1

    def get_lm_file_data(self, file_name):
        """Returns a dataframe containing the list mode data for only the specified file."""
        if file_name in self.lm_file_indices:
            indices = self.lm_file_indices[file_name]
            return self.lm_frame[indices[0]:indices[1]]
        else:
            raise ValueError(f'no file {file_name} in {self.name}.')

    def get_trace(self, time_id):
        """Returns trace data for the given time id."""
        if time_id in self.traces:
            return self.traces[time_id].copy(deep=True)
        else:
            raise ValueError(f"No trace with name '{time_id}' for scintillator '{self.name}'.")

    def get_trace_ids(self):
        """Returns a list of time ids for traces that are currently being stored."""
        return list(self.traces.keys())

    def find_matching_trace_id(self, count_time, first_sec, trace_list=None, count_file_index=None):
        """Returns the time id of the trace most likely to be associated with the given count."""
        # Checking to see that the count is inside the day or in the ~500 seconds of the next day included sometimes
        if count_time < 0 or count_time > params.SEC_PER_DAY + 500:
            return ''

        if trace_list is None:
            trace_list = self.get_trace_ids()

        # Getting the time range of the file the event occurred in
        if count_file_index is None:
            index = self.find_lm_file_index(count_time)
            if index != -1:
                count_file_extrema = self.lm_filetime_extrema[index]
            else:
                return ''

        else:
            count_file_extrema = self.lm_filetime_extrema[count_file_index]

        best_match = ''
        best_diff = float('inf')
        for time_id in trace_list:
            if time_id in self.traces:
                # No matter what, this time is probably wrong because it doesn't account
                # for any buffers except the first
                try:
                    # Datetime is sometimes a series object with two entries?
                    # Looks like a reader bug. Should always be a single datetime object
                    trace_time = self.traces[time_id]['DateTime'][0].timestamp() - first_sec
                except AttributeError:
                    continue

                # print(time_id, count_time, trace_time)
                if count_file_extrema[0] <= trace_time <= count_file_extrema[1]:
                    diff = abs(trace_time - count_time)
                    if diff < best_diff:
                        best_match = time_id
                        best_diff = diff

        return best_match
