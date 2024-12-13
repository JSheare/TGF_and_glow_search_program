"""A base class for keeping track of lightning data and associated information."""
import contextlib as contextlib
import datetime as dt
import gc as gc
import glob as glob
import json as json
import matplotlib.pyplot as plt
import numpy as np
import os as os
import pandas as pd
import psutil as psutil
import warnings

import tgfsearch.parameters as params
import tgfsearch.tools as tl
import tgfsearch.DataReaderTimetrack2 as Dr
from tgfsearch.detectors.scintillator import Scintillator


class Detector:
    """A class used to store all relevant information about an instrument and its data for a day.

    The Detector class is used to store the name of the detector, the date in various formats,
    and the actual data for the requested day in a single, centralized location.

    Parameters
    ----------
    unit : str
        The name of the instrument.
    date_str : str
        The timestamp for the requested day in yymmdd format.
    print_feedback : bool
        A flag specifying whether feedback should be printed to stdout or not.

    Attributes
    ----------
    log : _io.TextIO
        The file where actions and findings are logged.
    first_sec : float
        The first second of the day in EPOCH time.
    full_date_str : str
        The timestamp for the requested in day in yyyy-mm-dd format.
    dates_stored : list
        A list of dates currently being stored in the Detector.
    deployment : dict
        Deployment information for the instrument on the requested day (if available).
    default_scintillator : str
        A string representing the default scintillator.
    spectra_params : dict
        A dictionary containing various parameters used to make spectra for the Detector.
    default_data_loc : str
        The default directory for an instrument's raw data.
    _import_loc : str
        The directory where data files for the day are located.
    _results_loc : str
        The directory where results will be exported.
    _scintillators : dict
        A dictionary containing Scintillators. These keep track of data for each of the instrument's
        scintillators. Note the name mangling underscore.
    scint_list : list
        A list of the instrument's scintillator names.
    processed : bool
        A flag for whether the Detector should import processed data.

    """

    def __init__(self, unit, date_str, print_feedback=False):
        # Basic information
        self.unit = unit.upper()
        self.date_str = date_str  # yymmdd
        self.print_feedback = print_feedback
        self.log = None
        self.first_sec = tl.get_first_sec(self.date_str)
        self.full_date_str = dt.datetime.utcfromtimestamp(int(self.first_sec)).strftime('%Y-%m-%d')  # yyyy-mm-dd
        self.dates_stored = [date_str]
        self.deployment = self._get_deployment()
        self.default_scintillator = 'LP'  # Don't change this unless you have a really good reason

        # Detector-specific information
        self.spectra_params = {'bin_range': 0, 'bin_size': 0}
        self.default_data_loc = ''
        self._import_loc = ''
        self._results_loc = os.getcwd() + f'/Results/{self.unit}/{self.date_str}'
        self._scintillators = {self.default_scintillator: Scintillator(self.default_scintillator, '0')}
        self.scint_list = []

        self.processed = False

    def __str__(self):
        """String casting overload. Returns a string of the form 'Detector(unit, date_str)'."""
        return f'Detector({self.unit}, {self.date_str})'

    # Debugging string dunder
    def __repr__(self):
        """Debugging string dunder method. Returns a string of the form 'Detector(unit, date_str)' along
        with some info about which Scintillators have data."""
        scintillators_with_data = []
        has_data = False
        for scintillator in self._scintillators:
            if self._scintillators[scintillator]:
                has_data = True
                scintillators_with_data.append(scintillator)

        default_string = self.__str__()
        data_string = f' in {scintillators_with_data}' if has_data else ''
        return default_string + f' Has data = {has_data}' + data_string

    def __iter__(self):
        """Iterator dunder. Returns a generator that yields the Detector's scintillator names."""
        for scintillator in self.scint_list:
            yield scintillator

    def __bool__(self):
        """Bool casting overload. Returns True if data for the default scintillator is present."""
        return self.data_present_in(self.default_scintillator)

    def __add__(self, operand_detector):
        """Addition operator overload. Returns a new Detector containing data from the current Detector
        and the provided one. See splice() method documentation for more information."""
        return self.splice(operand_detector)

    def _get_deployment(self):
        """Returns a dictionary full of deployment information for the instrument on its specified date."""
        for file in glob.glob(
                f'{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/deployments/'
                f'{self.unit}_deployment_*.json'):
            file_dates = file.split('deployment')[-1][1:].split('.')[0].split('_')
            if int(file_dates[0]) <= int(self.date_str) <= int(file_dates[1]):
                with open(file, 'r') as deployment:
                    return json.load(deployment)

        return {'location': 'no location listed', 'instrument': self.unit, 'start_date': '000000', 'end_date': '000000',
                'utc_to_local': 0.0, 'dst_in_region': False, 'weather_station': '', 'sounding_station': '',
                'latitude': 0., 'longitude': 0., 'altitude': 0.,
                'notes': ''}

    def file_form(self, eRC):
        """Returns the regex for a scintillator's files given the scintillator's eRC serial number."""
        return ''

    def get_import_loc(self):
        """Returns the directory where data will be imported from.

        Returns
        -------
        str
            The import directory as a string.

        """

        return self._import_loc

    def set_import_loc(self, loc):
        """Sets the directory where data will be imported from.

        Parameters
        ----------
        loc : str
            The import directory as a string.

        """

        if isinstance(loc, str):
            if len(loc) > 0:
                if loc[-1] == '/':
                    loc = loc[:-1]

            self._import_loc = loc
        else:
            raise TypeError('loc must be a string.')

    def get_results_loc(self):
        """Returns the directory where all results will be stored.

        Returns
        -------
        str
            The results directory as a string.

        """

        return self._results_loc

    def set_results_loc(self, loc, subdir=True):
        """Sets the directory where all results will be stored.

        Parameters
        ----------
        loc : str
            The results directory as a string.

        subdir : bool
            Optional. If True, results will be exported to a subdirectory of the form 'Results/unit/yymmdd' inside
            the specified location rather than straight to the location itself. True by default.

        """

        if isinstance(loc, str):
            if len(loc) > 0:
                if loc[-1] == '/':
                    loc = loc[:-1]

            if subdir:
                self._results_loc = loc + f'/Results/{self.unit}/{self.date_str}'
            else:
                self._results_loc = loc
        else:
            raise TypeError('loc must be a string.')

    def use_processed(self, overwrite_import_loc=True):
        """Tells the Detector to import processed data instead of normal raw data. Only available for Godot."""
        if self.is_named('GODOT'):
            self.processed = True
            if overwrite_import_loc:
                self._import_loc = f'/media/godot/godot/monthly_processed/{self.date_str[0:4]}'

        else:
            raise ValueError('processed data mode is only available for GODOT.')

    def is_named(self, name):
        """Returns True if the Detector has the same name as the passed string.

        Parameters
        ----------
        name : str
            The name being queried.

        Returns
        -------
        bool
            True if the name matches the name of the Detector, False otherwise.

        """

        return True if name.upper() == self.unit else False

    def data_present_in(self, scintillator, data_type='lm'):
        """Returns True if data is present for the specified scintillator and False otherwise.

        Parameters
        ----------
        scintillator : str
            The scintillator's name. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        data_type : str
            Optional. The type of data to check for. Use 'lm' to check for list mode data and 'trace' to check
            for trace data. Checks for list mode data by default.

        Returns
        -------
        bool
            True if data is present in the specified scintillator, False otherwise.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].data_present(data_type)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def get_attribute(self, scintillator, attribute, deepcopy=True):
        """Returns the requested attribute for a particular scintillator.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        attribute : str
            The name of the attribute of interest.
        deepcopy : bool
            Optional. If True, a deepcopy of the requested attribute will be returned. True by default.

        Returns
        -------
        str || list || numpy.ndarray || dict || pandas.core.frame.DataFrame
            String if 'eRC' is requested; list if 'lm_filelist' or 'lm_file_ranges' is requested;
            numpy array  if 'time', 'energy', or 'wc' is requested; Reader if 'reader' is requested; dataframe
            if 'lm_frame' is requested, etc.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].get_attribute(attribute, deepcopy)

        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def set_attribute(self, scintillator, attribute, new_info, deepcopy=True):
        """Updates the requested attribute for a particular scintillator.
        Note: new info must be of the same type as the old.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        attribute : str
            The name of the attribute of interest.
        new_info : any
            The new information for the requested attribute.
        deepcopy : bool
            Optional. If True, the requested attribute will be set to a deepcopy of new_info. True by default.

        """

        if scintillator in self._scintillators:
            self._scintillators[scintillator].set_attribute(attribute, new_info, deepcopy)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def get_lm_data(self, scintillator, column, file_name=None):
        """Returns a single column of list mode data as a numpy array.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        column : str
            The column of interest.
        file_name : str
            Optional. The name of the file to get data for. If not specified,
            data will be retrieved for the whole day.

        Returns
        -------
        numpy.ndarray
            A numpy array containing the requested data column.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].get_lm_data(column, file_name)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def set_lm_data(self, scintillator, column, new_data, file_name=None):
        """Sets a single column of list mode data to the new data specified.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        column : str
            The column of interest.
        new_data : numpy.ndarray
            A numpy array containing the new data.
        file_name : str
            Optional. The name of the file to set data for. If not specified,
            data will be set for the whole day.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].set_lm_data(column, new_data, file_name)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def find_lm_file_index(self, scintillator, count_time):
        """Helper function. Returns the index of the list mode file that the given count occurred in.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator that the count is from. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        count_time : float
            The time that the count occurred at (in seconds of day).

        Returns
        -------
        int
            The index of the list mode file that the given count occurred in. Returns -1 if the count
            isn't in any of the files.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].find_lm_file_index(count_time)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def find_lm_file(self, scintillator, count_time):
        """Returns the name of the list mode file that the given count occurred in.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator that the count is from. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        count_time : float
            The time that the count occurred at (in seconds of day).

        Returns
        -------
        str
            The name of the list mode file that the given count occurred in. Returns an empty string if the count
            isn't in any of the files.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].find_lm_file(count_time)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def get_lm_file(self, scintillator, file_name, deepcopy=True):
        """Returns the list mode data for the specified list mode file.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator that the file is from. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        file_name : str
            The name of the file that data is being requested for. Note that this must be the *full* name of the file,
            including the path from the root directory.
        deepcopy : bool
            Optional. If True, a deepcopy of the file frame will be returned. True by default.

        Returns
        -------
        pandas.core.frame.DataFrame
            A dataframe with the data for the requested list mode file.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].get_lm_file(file_name, deepcopy)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def get_trace(self, scintillator, trace_name, deepcopy=True):
        """Returns the trace data for the given scintillator and trace name.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        trace_name : str
            The name of the trace file.
        deepcopy : bool
            Optional. If True, a deepcopy of the trace frame will be returned. True by default.

        Returns
        -------
        pandas.core.frame.DataFrame
            A dataframe containing the trace data for the given name.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].get_trace(trace_name, deepcopy)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def get_trace_names(self, scintillator):
        """Returns a list of names of the traces that are currently being stored.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.

        Returns
        -------
        list
            A list of names of the traces that are currently being stored.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].get_trace_names()
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def find_matching_traces(self, scintillator, count_time, trace_list=None):
        """Finds the traces that could be a match for the given count (if they exist).

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        count_time : float
            The time that the count occurred at (in seconds of day).
        trace_list : list
            Optional. A list of traces that could be valid matches. If not provided, all traces stored for the
            scintillator of interest will be considered valid.

        Returns
        -------
        list
            A list of trace names that could be matches.

        """

        if scintillator in self._scintillators:
            return self._scintillators[scintillator].find_matching_traces(count_time, self.date_str, trace_list)
        else:
            raise ValueError(f"'{scintillator}' is not a valid scintillator.")

    def get_fileset_size(self):
        """Returns the total size (in bytes) of all the currently listed files for the day."""
        total_file_size = 0
        for scintillator in self._scintillators:
            for file in self.get_attribute(scintillator, 'lm_filelist', deepcopy=False):
                total_file_size += tl.file_size(file)

            for file in self.get_attribute(scintillator, 'trace_filelist', deepcopy=False):
                total_file_size += tl.file_size(file)

        return total_file_size

    def clear(self, clear_filelists=True):
        """Clears all data currently stored in the Detector.

        Parameters
        ----------
        clear_filelists : bool
            Optional. If True, lists of files stored in the detector will be also be cleared. True by default.

        """

        for scintillator in self._scintillators:
            self._scintillators[scintillator].clear(clear_filelists)

        if clear_filelists:
            self.dates_stored = [self.date_str]

    def _import_lm_data(self, scintillator):
        """Imports list mode data for the given scintillator."""
        lm_filelist = self.get_attribute(scintillator, 'lm_filelist')
        if len(lm_filelist) < 1:
            if self.log is not None:
                print('Missing list mode data', file=self.log)
                print('', file=self.log)

            if self.print_feedback:
                print('Missing list mode data')

            return

        file_frames = []
        file_ranges = []
        file_indices = {}
        file_time_gaps = []
        prev_second = 0
        start_index = 0
        files_imported = 0

        if self.log is not None:
            print('List Mode Files:', file=self.log)
            print('File|Import Success|File Time Gap (sec)', file=self.log)

        if self.print_feedback:
            print(f'Importing {len(lm_filelist)} list mode files...')

        # Importing the data
        for file in lm_filelist:
            if self.processed:
                energy, time = np.loadtxt(file, skiprows=1, usecols=(0, 2), unpack=True)
                data = pd.DataFrame.from_dict({'time': time, 'energy': energy})
                file_frames.append(data)
            else:
                reader = self.get_attribute(scintillator, 'reader', deepcopy=False)
                # Try-except block to handle reader errors
                try:
                    # The first with disables prints from the data reader; The second with suppresses annoying
                    # numpy warnings
                    with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore', category=RuntimeWarning)
                            data = reader.read(file)

                except Exception as ex:
                    # Files that generate reader errors are skipped
                    if self.log is not None:
                        print(f'{file}|False|N/A', file=self.log)
                        print(f'    Error importing file: {ex}', file=self.log)

                    continue

                if 'energies' in data.columns:
                    data.rename(columns={'energies': 'energy'}, inplace=True)

            # first_second = data['SecondsOfDay'].iloc[0]
            # last_second = data['SecondsOfDay'].iloc[-1]

            # Above would be better, but counts are very occasionally out of chronological order for some reason
            first_second = data['SecondsOfDay'].min()
            last_second = data['SecondsOfDay'].max()

            # Determines the time gaps between adjacent files
            file_time_gap = first_second - prev_second if files_imported > 0 else 0.0
            file_time_gaps.append(file_time_gap)
            prev_second = last_second
            file_ranges.append([first_second, last_second])
            file_frames.append(data)
            if self.log is not None:
                print(f'{file}|True|{file_time_gap}', file=self.log)

            # Keeps track of file indices in the larger dataframe
            data_length = len(data.index)
            file_indices[file] = [start_index, start_index + data_length]
            start_index += data_length

            files_imported += 1

        if self.print_feedback:
            print(f'{files_imported}/{len(lm_filelist)} list mode files imported')

        if len(file_frames) > 0:
            # Correcting for the fact that the first 200-300 seconds of the next day are usually included
            # in the last file
            last_times = file_frames[-1]['SecondsOfDay'].to_numpy()
            # Times belonging to the next day will always be < 500
            day_change = np.where(np.diff(last_times) < -(params.SEC_PER_DAY - 500))[0]
            if len(day_change) > 0:
                change_index = int(day_change[0]) + 1
                for i in range(change_index, len(last_times)):
                    last_times[i] += params.SEC_PER_DAY

                file_frames[-1]['SecondsOfDay'] = last_times

            # Correcting the last file's time ranges too
            if (file_ranges[-1][1] - file_ranges[-1][0]) < -(params.SEC_PER_DAY - 500):
                file_ranges[-1][1] += params.SEC_PER_DAY

            # Makes the final dataframe and stores it
            all_data = pd.concat(file_frames, axis=0)

            self.set_attribute(scintillator, 'lm_frame', all_data, deepcopy=False)
            self.set_attribute(scintillator, 'lm_file_ranges', file_ranges, deepcopy=False)
            self.set_attribute(scintillator, 'lm_file_indices', file_indices, deepcopy=False)

            if self.log is not None:
                print('', file=self.log)
                print(f'Total Counts: {len(all_data.index)}', file=self.log)
                print(f'Average time gap: {sum(file_time_gaps) / len(file_time_gaps)}', file=self.log)
                print('', file=self.log)
        else:
            if self.log is not None:
                print('', file=self.log)

    def _import_trace_data(self, scintillator):
        """Imports trace data for the given scintillator."""
        trace_filelist = self.get_attribute(scintillator, 'trace_filelist')
        if len(trace_filelist) < 1:
            if self.log is not None:
                print('No trace data', file=self.log)
                print('', file=self.log)

            if self.print_feedback:
                print('No trace data')

            return

        traces = {}
        files_imported = 0

        if self.log is not None:
            print('Trace Files:', file=self.log)
            print('File|Import Success|', file=self.log)

        if self.print_feedback:
            print(f'Importing {len(trace_filelist)} trace files...')

        # Importing the data
        for file in trace_filelist:

            # Try-except block for handling reader errors
            try:
                # The first with disables prints from the data reader; The second with suppresses annoying
                # numpy warnings
                with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', category=RuntimeWarning)
                        data = Dr.fileNameToData(file, {})

            except Exception as ex:
                # Files that generate reader errors are skipped
                if self.log is not None:
                    print(f'{file}|False|', file=self.log)
                    print(f'    Error importing file: {ex}', file=self.log)

                continue

            traces[file] = data
            if self.log is not None:
                print(f'{file}|True|', file=self.log)

            files_imported += 1

        if self.print_feedback:
            print(f'{files_imported}/{len(trace_filelist)} trace files imported')

        # Storing the traces
        if len(traces) > 0:
            self.set_attribute(scintillator, 'traces', traces, deepcopy=False)

        if self.log is not None:
            print('', file=self.log)

    def import_data(self, existing_filelists=False, import_traces=True, import_lm=True, mem_frac=1.):
        """Imports data from data files into arrays and then updates them into the detector's
        scintillator objects.

        Parameters
        ----------
        existing_filelists : bool
            Optional. If True, the function will use the file lists already stored in the Detector.
        import_traces : bool
            Optional. If True, the function will import any trace files it finds. True by default.
        import_lm : bool
            Optional. If True, the function will import any list mode data it finds. True by default.
        mem_frac : float
            Optional: The maximum fraction of currently available system memory that the Detector is allowed to use for
             data (not including overhead). If the dataset is larger than this limit, a MemoryError will be raised.
             1.0 (100% of system memory) by default.

        """

        if len(self.dates_stored) > 1:
            raise RuntimeError("cannot import multiple days' data.")

        if not existing_filelists:
            # Locates the files to be imported
            for scintillator in self._scintillators:
                eRC = self.get_attribute(scintillator, 'eRC')
                # Here in case the data files are grouped into daily folders
                try:
                    complete_filelist = glob.glob(f'{self._import_loc}/{self.file_form(eRC)}'.replace('\\', '/'))
                    assert len(complete_filelist) > 0, 'Empty filelist'

                except AssertionError:
                    complete_filelist = glob.glob(f'{self._import_loc}/{self.date_str}'
                                                  f'/{self.file_form(eRC)}'.replace('\\', '/'))

                lm_filelist, trace_filelist = tl.separate_data_files(tl.filter_data_files(complete_filelist))
                if import_lm:
                    self._scintillators[scintillator].lm_filelist = lm_filelist

                if import_traces:
                    self._scintillators[scintillator].trace_filelist = trace_filelist

        # Checking to make sure that currently listed data won't go over the memory limit
        if self.get_fileset_size() > psutil.virtual_memory()[1] * mem_frac:
            raise MemoryError('dataset larger than specified limit.')

        if self.log is not None:
            print('', file=self.log)

        for scintillator in self._scintillators:
            eRC = self.get_attribute(scintillator, 'eRC')
            if self.log is not None:
                print(f'For eRC {eRC} ({scintillator}):', file=self.log)

            if self.print_feedback:
                print('')
                print(f'For eRC {eRC} ({scintillator}):')

            # Importing list mode data
            if import_lm:
                self._import_lm_data(scintillator)

            # Importing trace data
            if import_traces:
                self._import_trace_data(scintillator)

        gc.collect()

    def make_spectra(self, scintillator):
        """Makes an energy spectra histogram for the requested scintillator.

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            Two numpy arrays: the energy bins, and the energy histogram.

        """

        if self.data_present_in(scintillator):
            energy_bins = np.arange(0.0, self.spectra_params['bin_range'], self.spectra_params['bin_size'])
            data = self.get_attribute(scintillator, 'lm_frame', deepcopy=False)
            energy_hist = np.histogram(data['energy'], bins=energy_bins)[0]
            return energy_bins[:-1], energy_hist  # Bins is always longer than hist by one
        else:
            raise ValueError(f"data for '{scintillator}' is either missing or hasn't been imported yet.")

    def plot_spectra(self, scintillator):
        """Plots and saves the energy spectra histogram for the given scintillator

        Parameters
        ----------
        scintillator : str
            The name of the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.

        """

        figure = plt.figure(figsize=[20, 11.0])
        ax = figure.add_subplot()
        figure.suptitle(f'Energy Spectrum for {scintillator}, {self.full_date_str}')
        ax.set_xlabel('Energy Channel')
        ax.set_ylabel('Counts / bin')
        ax.set_yscale('log')
        energy_bins, energy_hist = self.make_spectra(scintillator)
        ax.bar(energy_bins, energy_hist, color='r',
               width=self.spectra_params['bin_size'] / 2, zorder=1)

        # Saves the figure
        figure.savefig(f'{self._results_loc}/{scintillator}_Spectrum.png', dpi=500)
        figure.clf()
        plt.close(figure)
        gc.collect()

    def splice(self, operand_detector):
        """Returns a new Detector with the combined data of the current Detector and the one provided.

        Parameters
        ----------
        operand_detector : tgfsearch.detectors.detector.Detector
            Another Detector.

        Returns
        -------
        tgfsearch.detectors.detector.Detector
            A new Detector containing the data of the current Detector and the one provided.

            Three things to note:
                - This Detector will store the date of the *earlier* operand Detector

                - All time-related data will be adjusted to reflect the difference in date between the operand Detectors

                - Using the method import() with this Detector won't work. Trying this will result in a RuntimeError

        """

        if operand_detector.unit == self.unit:
            if int(self.date_str) < int(operand_detector.date_str):
                new_detector = type(self)(self.unit, self.date_str, print_feedback=self.print_feedback)
                earlier = self
                later = operand_detector
            elif int(self.date_str) > int(operand_detector.date_str):
                new_detector = type(self)(operand_detector.unit, operand_detector.date_str,
                                          print_feedback=operand_detector.print_feedback)
                earlier = operand_detector
                later = self
            else:
                raise ValueError('cannot splice the same day into itself.')

            new_detector.dates_stored.append(later.date_str)

            # Measuring the number of days between the earlier and later date
            day_difference = 0
            rolled_date = earlier.date_str
            while rolled_date != later.date_str:
                rolled_date = tl.roll_date_forward(rolled_date)
                day_difference += 1

            for scintillator in self._scintillators:
                # Combining list mode file lists
                new_detector.set_attribute(scintillator, 'lm_filelist',
                                           earlier.get_attribute(scintillator, 'lm_filelist') +
                                           later.get_attribute(scintillator, 'lm_filelist'),
                                           deepcopy=False)

                # Combining list mode file ranges
                # Updating ranges from later Detector to reflect the date difference
                new_ranges = later.get_attribute(scintillator, 'lm_file_ranges')
                for new_range in new_ranges:
                    new_range[0] += day_difference * params.SEC_PER_DAY
                    new_range[1] += day_difference * params.SEC_PER_DAY

                new_detector.set_attribute(scintillator, 'lm_file_ranges',
                                           earlier.get_attribute(scintillator, 'lm_file_ranges') +
                                           new_ranges, deepcopy=False)

                # Combining list mode file indices
                # Updating indices from later Detector to reflect their new positions in the data frame
                new_start = len(earlier.get_attribute(scintillator, 'lm_frame', deepcopy=False))
                new_indices = later.get_attribute(scintillator, 'lm_file_indices')
                for file in new_indices:
                    new_indices[file][0] += new_start
                    new_indices[file][1] += new_start

                new_indices.update(earlier.get_attribute(scintillator, 'lm_file_indices'))
                new_detector.set_attribute(scintillator, 'lm_file_indices', new_indices, deepcopy=False)

                # Combining list mode data frames
                # Updating later frame to reflect the difference in days
                later_frame = later.get_attribute(scintillator, 'lm_frame')
                later_frame['SecondsOfDay'] += day_difference * params.SEC_PER_DAY

                new_lm_frame = pd.concat([
                    earlier.get_attribute(scintillator, 'lm_frame', deepcopy=False),
                    later_frame], axis=0)

                new_detector.set_attribute(scintillator, 'lm_frame', new_lm_frame, deepcopy=False)

                # Combining trace file lists
                new_detector.set_attribute(scintillator, 'trace_filelist',
                                           earlier.get_attribute(scintillator, 'trace_filelist') +
                                           later.get_attribute(scintillator, 'trace_filelist'),
                                           deepcopy=False)

                # Combining trace tables
                new_traces = earlier.get_attribute(scintillator, 'traces')
                new_traces.update(later.get_attribute(scintillator, 'traces'))
                new_detector.set_attribute(scintillator, 'traces', new_traces, deepcopy=False)

            return new_detector
        else:
            raise TypeError(f"cannot splice '{self.unit}' with '{operand_detector.unit}'.")
