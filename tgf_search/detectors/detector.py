"""A base class for keeping track of data and associated information."""
import glob as glob
import os as os
import contextlib as contextlib
import psutil as psutil
import numpy as np
import pandas as pd
import json as json
import scipy.signal as signal
import matplotlib.pyplot as plt
import datetime as dt
from matplotlib.widgets import Slider

import tgf_search.DataReaderFinal as Dr
# import tgf_search.DataReaderTimetrack2 as Dr
import tgf_search.tools as tl
import tgf_search.parameters as params
from tgf_search.detectors.scintillator import Scintillator


class Detector:
    """Used to store all relevant information about the detector and the data.

    The detector object is used to store the name of the detector, the requested date for analysis in various formats,
    and the actual data in a single, centralized location.

    Parameters
    ----------
    unit : str
        The name of the detector that the analysis is being requested for.
    first_sec : float
        The first second in EPOCH time of the day that data analysis is being requested for.
    modes : list
        A list of the requested modes that the program should operate under.
    print_feedback : bool
        A flag specifying whether the object should print feedback to stdout or not.

    Attributes
    ----------
    log : file
        The .txt file where actions and findings are logged.
    date_str : str
        The timestamp for the requested day in yymmdd format.
    full_date_str : str
        The timestamp for the requested in day in yyyy-mm-dd format.
    location : dict
        Location information for the detector on the requested day.
    lp_calibrated : bool
        A flag for whether the large plastic scintillator has been calibrated.
    default_scintillator : str
        A string representing the default scintillator. Data from this scintillator must be present for the program
        to run.
    long_event_scint_list : list
        A list of the scintillators to used by the search script's long event search algorithm for a detector.
    calibration_params : dict
        A dictionary containing various parameters used in the detector calibration algorithm.
    default_data_loc : str
        The default directory for a detector's raw data.
    import_path : str
        The directory where the requested data files are located.
    results_loc : str
        The directory where program results will be exported.
    regex : function
        A lambda function that, when given the eRC serial number, returns the regex for a scintillator's files.
    scintillators : dict
        A dictionary containing scintillator objects. These objects keep track of data for each scintillator.
    scint_list : list
        A list of the detector's scintillators.
    processed : bool
        A flag for whether the program should operate in "processed" mode or not. Under this mode, the program will use
        the built-in file path for GODOT processed data on Sol.
    template : bool
        A flag for whether the program should operate in "template" mode or not. This mode allows the user to generate
        and tweak templates that the program uses to calibrate the large plastic scintillator. These must be made for
        each new detector location.
    gui : bool
        A flag for whether the object is being instantiated from the GUI (usually via the search script).

    """

    def __init__(self, unit, first_sec, modes, print_feedback=False):
        # Basic information
        self.unit = unit
        self.first_sec = first_sec
        self.modes = modes
        self.print_feedback = print_feedback
        self.log = None
        self.date_str = dt.datetime.utcfromtimestamp(int(first_sec)).strftime('%y%m%d')  # In format yymmdd
        self.full_date_str = dt.datetime.utcfromtimestamp(int(first_sec)).strftime('%Y-%m-%d')  # In format yyyy-mm-dd
        self.location = None
        self.lp_calibrated = False
        self.default_scintillator = 'LP'  # Don't change this unless you have a really good reason

        # Specific detector information
        self.long_event_scint_list = []
        self.calibration_params = {'bin_range': 0, 'bin_size': 0, 'band_starts': [0, 0],
                                   'band_ends': [0, 0], 'template_bin_plot_edge': 0}

        self.default_data_loc = ''
        self.import_path = ''
        self.results_loc = os.getcwd() + '/'
        self.regex = lambda eRC: ''
        self.scintillators = {'NaI': Scintillator('NaI', '0'), 'LP': Scintillator('LP', '0')}

        self.scint_list = []

        # Modes
        self.processed = False
        self.template = False
        self.gui = False

    # String casting overload
    def __str__(self):
        return f'Detector({self.unit}, {self.first_sec}, {self.modes})'

    # Debugging string dunder
    def __repr__(self):
        scintillators_with_data = []
        has_data = False
        for scintillator in self:
            if self.scintillators[scintillator]:
                has_data = True
                scintillators_with_data.append(scintillator)
                break

        default_string = self.__str__()
        data_string = f' in {scintillators_with_data}' if has_data else ''
        return default_string + f' Has data = {has_data}' + data_string

    # Iterator dunder (for use in loops)
    def __iter__(self):
        for scintillator in self.scint_list:
            yield scintillator

    # Bool casting overload. Returns True if data for default scintillator is present
    def __bool__(self):
        return self.is_data_present(self.default_scintillator)

    def check_processed(self):
        """Checks to see if processed is one of the user specified modes, and raises an error if the detector isn't
        Godot."""
        if 'processed' in self.modes:
            if self.is_named('GODOT'):
                self.processed = True
                self.import_path = f'/media/godot/godot/monthly_processed/{self.date_str[0:4]}'
            else:
                raise ValueError('ValueError: processed data mode is only available for GODOT.')

    def check_gui(self):
        """Checks to see if the object is being instantiated from the GUI (via the search script usually), and
        changes the import and export directories if the user specified different ones from the default."""
        self.gui = True if 'GUI' in self.modes else False
        if self.gui:
            if self.modes[-2] != 'none':
                if self.modes[-2] != '/':
                    self.results_loc = self.modes[-2] + '/'
                else:
                    self.results_loc = self.modes[-2]

            if self.modes[-1] != 'none':
                if self.modes[-1] == '/':
                    self.import_path = self.modes[-1][:-1]
                else:
                    self.import_path = self.modes[-1]

    def get_location(self, deployment_file_loc):
        """Returns a dictionary full of location information for a detector on it's specified date."""
        for file in glob.glob(f'{deployment_file_loc}/{self.unit.lower()}_*_*.json'):
            if int(file[6:12]) <= int(self.date_str) <= int(file[13:19]):
                with open(file, 'r') as deployment:
                    return json.load(deployment)

        return {'Location': 'no location listed', 'Instrument': self.unit, 'Start date': '', 'End date': '',
                'UTC conversion to local time': '', 'Nearest weather station': '', 'Daylight Savings?': '',
                'Latitude (N)': '', 'Longitude (E, 0-360)': '', 'Altitude (km)': '',
                'Notes': ''}

    def is_named(self, name):
        """Returns True if the detector has the same name as the passed string.

        Parameters
        ----------
        name : str
            The potential detector name being queried.

        Returns
        -------
        bool
            True if the name matches the name of the detector, False otherwise.
        """
        return True if name.upper() == self.unit else False

    def is_data_present(self, scintillator):
        """Returns True if data is present for the requested scintillator and False otherwise.

        Parameters
        ----------
        scintillator : str
            A string corresponding to the scintillator of interest.

        Returns
        -------
        bool
            True if data present in the requested scintillator, False otherwise.
        """

        return bool(self.scintillators[scintillator])

    def get_attribute(self, scintillator, attribute):
        """Retrieves the requested attribute for a particular scintillator

        \n
        Attribute summary:
            eRC: the scintillator's serial number.

            filelist: the list of files for the day.

            filetime_extrema: a list of lists. Each sublist contains the first and last second present in every file.

            calibration: a list of two energies (in Volts) used to calibrate the scintillator.

            time: a numpy array of each count's time (in seconds since start of day).

            energy: a numpy array of each count's energy (in Volts).

            wc: a numpy array of each count's wallclock time.

            passtime: timing information for the previous file imported.

        Parameters
        ----------
        scintillator : str
            A string corresponding to the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        attribute : str
            A string corresponding to the scintillator attribute of interest. See attribute
            summary above for a full list of allowed values.

        Returns
        -------
        str, list, np.array, dict
            String if 'eRC' is requested; list if 'filelist', 'calibration', or 'filetime_extrema' is requested;
            numpy array  if 'time', 'energy', or 'wc' is requested; dictionary if 'passtime' is requested.

        """

        if scintillator in self.scintillators:
            medium = self.scintillators[scintillator]
            return medium.get_attribute(attribute)

        else:
            raise ValueError('ValueError: not a valid scintillator.')

    def update_attribute(self, scintillator, attribute, new_info):
        """Updates the requested attribute for a particular scintillator.

        \n
        Attribute summary:
            eRC: the scintillator's serial number.

            filelist: the list of files for the day.

            filetime_extrema: a list of lists. Each sublist contains the first and last second present in every file.

            calibration: a list of two energies (in Volts) used to calibrate the scintillator.

            time: a numpy array of each count's time (in seconds since start of day).

            energy: a numpy array of each count's energy (in Volts).

            wc: a numpy array of each count's wallclock time.

            passtime: timing information for the previous file imported.

        Parameters
        ----------
        scintillator : str
            A string corresponding to the scintillator of interest. Allowed values (detector dependent):
            'NaI', 'SP', 'MP', 'LP'.
        attribute : str/list
            For only one attribute: a string corresponding to the scintillator attribute of interest. See attribute
            summary above for a full list of allowed values.
            For multiple attributes: a list of strings corresponding to the scintillator attributes of interest.
        new_info : any/list
            For only one attribute: the new information for the requested attribute.
            For multiple attributes: a list of the new information for each requested attribute.

        """

        if scintillator in self:
            medium = self.scintillators[scintillator]
            if type(attribute) is list and type(new_info) is list:
                if len(attribute) != len(new_info):
                    raise ValueError('ValueError: attribute and new info must be the same length.')

                for i in range(len(attribute)):
                    medium.update_attribute(attribute[i], new_info[i])

            elif type(attribute) is str:
                medium.update_attribute(attribute, new_info)

            else:
                raise TypeError('TypeError: if attribute is a list, new_info must also be a list.')

        else:
            raise ValueError('ValueError: not a valid scintillator.')

    def calculate_fileset_size(self):
        """Returns the total size (in bytes) of all the currently held files for the day."""
        total_file_size = 0
        for scintillator in self:
            filelist = self.get_attribute(scintillator, 'filelist')
            for file in filelist:
                total_file_size += os.path.getsize(file)

        return total_file_size

    def return_passtime(self):
        """Returns a dictionary containing the passtime attributes for each of the detector's scintillators.

        Returns
        ----------
        dict
            A dictionary containing the passtime dictionaries for each of the detector's scintillators.

        """

        passtime_dict = {}
        for scintillator in self:
            passtime = self.get_attribute(scintillator, 'passtime')
            passtime_dict.update({scintillator: passtime})

        return passtime_dict

    def update_passtime(self, passtime_dict):
        """Updates the passtime attributes for each of the detector's scintillators.

        Parameters
        ----------
        passtime_dict : dict
            A dictionary containing the passtime dictionaries for each of the detector's scintillators.

        """

        for scintillator in self:
            self.update_attribute(scintillator, 'passtime', passtime_dict[scintillator])

    def import_data(self, existing_filelists=False):
        """Imports data from data files into arrays and then updates them into the detector's
        scintillator objects.

        Parameters
        ----------
        existing_filelists : bool
            Optional. If True, the function will use the file lists already stored in the detector's scintillator
            objects.

        """

        for scintillator in self.scintillators:
            if existing_filelists:
                break

            eRC = self.get_attribute(scintillator, 'eRC')
            # Here in case the data files in a custom location are grouped into daily folders
            try:
                complete_filelist = glob.glob(f'{self.import_path}/{self.regex(eRC)}')
                assert len(complete_filelist) > 0, 'Empty filelist'

            except AssertionError:
                complete_filelist = glob.glob(f'{self.import_path}/{self.date_str}'
                                              f'/{self.regex(eRC)}')

            filelist = tl.filter_files(complete_filelist)
            self.update_attribute(scintillator, 'filelist', filelist)

        # Checks to see if the necessary files for a full search are present
        if len(self.get_attribute(self.default_scintillator, 'filelist')) == 0:
            missing_data_scints = []
            for scintillator in self.scintillators:
                if len(self.get_attribute(scintillator, 'filelist')) == 0:
                    missing_data_scints.append(scintillator)

            print('\n', file=self.log)
            print('No/missing necessary data for specified day.', file=self.log)
            print(f'Data missing in the following: {", ".join(missing_data_scints)}', file=self.log)
            if self.print_feedback:
                print('\n')
                print('No/missing necessary data for specified day.')
                print('\n')

            raise FileNotFoundError

        # Determines whether there is enough free memory to load the entire dataset
        total_file_size = self.calculate_fileset_size()
        available_memory = psutil.virtual_memory()[1] * params.TOTAL_MEMORY_ALLOWANCE_FRAC
        if (params.OPERATING_MEMORY_ALLOWANCE + total_file_size) > available_memory:
            raise MemoryError('MemoryError: not enough free memory to hold complete dataset.')

        for scintillator in self.scintillators:
            eRC = self.get_attribute(scintillator, 'eRC')
            filelist = self.get_attribute(scintillator, 'filelist')
            print('\n', file=self.log)
            print(f'For eRC {eRC} ({scintillator}):', file=self.log)
            if self.print_feedback:
                print('\n')
                print(f'For eRC {eRC} ({scintillator}):')

            try:
                # Tests to make sure that filelist isn't empty
                assert len(filelist) > 0, 'No files for this scintillator today'

                # Starts actually importing the data
                energy_list = []
                time_list = []
                wallclock_list = []
                filetime_extrema_list = []

                filetimes = np.array([])
                file_time_gaps = np.array([])
                last_second = 0.0
                files_imported = 0

                print('File|File Behavior|File Time Gap (sec)', file=self.log)
                filecount_switch = True
                for file in filelist:
                    if not self.gui and self.print_feedback:
                        print(f'{files_imported}/{len(filelist)} files imported', end='\r')
                    elif self.gui and filecount_switch and self.print_feedback:
                        print(f'Importing {len(filelist)} files...')
                        filecount_switch = False

                    # Try-except block to log files where GPS and wallclock disagree significantly
                    file_behavior = 'Normal'
                    try:
                        if self.processed:
                            e, t = np.loadtxt(file, skiprows=1, usecols=(0, 2), unpack=True)
                            energy_list.append(e)
                            time_list.append(t)
                            filetimes = t
                        else:
                            with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):  # disables prints from dr
                                # data, passtime = Dr.fileNameToData(file,
                                #                                    self.get_attribute(scintillator, 'passtime'))
                                data = Dr.fileNameToData(file)

                            # self.update_attribute(scintillator, 'passtime', passtime)
                            if 'energies' in data.columns:
                                energy_list.append(data['energies'].to_numpy())
                            else:
                                energy_list.append(data['energy'].to_numpy())

                            time_list.append(data['SecondsOfDay'].to_numpy())
                            wallclock_list.append(data['wc'].to_numpy())
                            filetimes = data['SecondsOfDay'].to_numpy()
                            filetime_extrema_list.append([filetimes[0], filetimes[-1]])

                    except Exception as ex:
                        if str(ex) == 'wallclock and GPS clocks in significant disagreement':
                            file_behavior = 'Disagreement'
                            pass
                        else:  # Mostly here so that if the reader ever runs into other errors I'll know about them
                            raise Exception('Reader Error')

                    # Determines the time gaps between adjacent files
                    first_second = filetimes[0]
                    file_time_gap = first_second - last_second if files_imported > 0 else 0
                    file_time_gaps = np.append(file_time_gaps, file_time_gap)
                    last_second = filetimes[-1]
                    files_imported += 1

                    print(f'{file}|{file_behavior}|{file_time_gap}', file=self.log)

                if self.print_feedback:
                    print(f'{files_imported}/{len(filelist)} files imported', end='\r')

                # Makes the final arrays and exports them
                times = np.concatenate(time_list)
                # Corrects for the fact that the first 200-300 seconds of the next day are included in the last file
                day_change_array = np.array(np.where(np.diff(times) < -80000))
                if day_change_array.size > 0:
                    change_index = int(day_change_array[0]) + 1
                    times = np.append(times[0:change_index], times[change_index:] + 86400.0)

                # Does it for the file time extrema too
                for k in range(int(len(filetime_extrema_list) / 8)):  # Last eighth of the files
                    last_file_extrema = filetime_extrema_list[-(k + 1)]
                    for j in range(2):
                        extrema = last_file_extrema[j]
                        if extrema < 500:  # Extrema belonging to the next day will always be < 500
                            last_file_extrema[j] = extrema + 86400

                    filetime_extrema_list[-(k + 1)] = last_file_extrema

                updated_attributes = ['time', 'energy', 'wc', 'filetime_extrema']
                updated_info = [times, np.concatenate(energy_list), np.concatenate(wallclock_list),
                                filetime_extrema_list]
                self.update_attribute(scintillator, updated_attributes, updated_info)

                print('\n', file=self.log)
                print(f'Total Counts: {len(np.concatenate(time_list))}', file=self.log)
                print(f'Average time gap: {np.sum(file_time_gaps) / len(file_time_gaps)}', file=self.log)
                print('\n', file=self.log)

            except AssertionError:
                print('Missing data for the specified day.', file=self.log)
                if self.print_feedback:
                    print('Missing data for the specified day.', end='\r')

                continue

            except Exception as ex:
                if str(ex) == 'Reader Error':
                    print('Error with data reader.', file=self.log)
                    if self.print_feedback:
                        print('Error with data reader.', end='\r')

                    continue
                else:
                    raise

    def make_spectra_hist(self, existing_spectra_dict):
        """Makes the energy spectra histograms for a chunk of the day (no calibration). Returns the histograms in a
        dictionary.

        Parameters
        ----------
        existing_spectra_dict : dict
            A dictionary containing energy spectra histograms for previous chunks of the day.

        Returns
        -------
        dict
            An updated version of existing_spectra_dict featuring the current chunk's contribution to the energy
            spectra histograms.

        """

        bin_range = self.calibration_params['bin_range']
        bin_size = self.calibration_params['bin_size']

        energy_bins = np.arange(0.0, bin_range, bin_size)
        for scintillator in self:
            energies = self.get_attribute(scintillator, 'energy')
            chunk_hist, bin_edges = np.histogram(energies, bins=energy_bins)
            if len(existing_spectra_dict[scintillator]) == 0:
                existing_spectra_dict[scintillator] = chunk_hist
            else:
                existing_spectra_dict[scintillator] = existing_spectra_dict[scintillator] + chunk_hist

        return existing_spectra_dict

    def _generate_hist(self, energy_bins, scintillator, existing_spectra=None):
        """Returns the energy spectra histogram for the requested scintillator."""
        energies = self.get_attribute(scintillator, 'energy')
        if existing_spectra is not None:
            energy_hist = existing_spectra[scintillator]
        else:
            energy_hist, bin_edges = np.histogram(energies, bins=energy_bins)

        return energy_hist

    def _make_template(self, energy_bins, energy_hist):
        """Makes a template that can be used in the LP calibration algorithm's cross-correlation."""
        if self.location['Location'] == 'no location listed':
            return

        bin_plot_edge = len(energy_bins) - 1  # Histogram array is shorter than bin array by 1 (no idea why)

        template_bin_plot_edge = self.calibration_params['template_bin_plot_edge']

        print('Entering template mode...', file=self.log)
        if self.print_feedback:
            print('Entering template mode...')
            print('\n')
            print('Use the sliders to adjust the line positions. The K40 line comes first.')

        def line_locs(e1, e2):
            return energy_bins[np.array([e1, e2]).astype(int)]

        # Initial vertical line positions
        edge1 = 0
        edge2 = 0

        fig, ax = plt.subplots()
        ax.set_xlabel('Energy Channel')
        ax.set_ylabel('Counts/bin')
        ax.set_yscale('log')
        ax.bar(energy_bins[0:template_bin_plot_edge], energy_hist[0:template_bin_plot_edge], color='r',
               width=self.calibration_params['bin_size'] / 2, zorder=1)
        lines = ax.vlines(line_locs(edge1, edge2), 0, np.amax(energy_hist), zorder=2, alpha=0.75)

        fig.subplots_adjust(bottom=0.30)

        # Slider for the Potassium 40 line
        ax1 = fig.add_axes([0.25, 0.15, 0.65, 0.03])
        edge1_slider = Slider(
            ax=ax1,
            label='K40',
            valmin=0,
            valmax=len(energy_bins[0:template_bin_plot_edge]) - 1,
            valinit=edge1,
            valstep=1,
        )

        # Slider for the Thorium line
        ax2 = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        edge2_slider = Slider(
            ax=ax2,
            label='T',
            valmin=0,
            valmax=len(energy_bins[0:template_bin_plot_edge]) - 1,
            valinit=edge2,
            valstep=1,
        )

        def update(val):
            nonlocal lines
            lines.remove()
            lines = ax.vlines(line_locs(edge1_slider.val, edge2_slider.val), 0, np.amax(energy_hist),
                              zorder=2, alpha=0.75)
            fig.canvas.draw_idle()

        edge1_slider.on_changed(update)
        edge2_slider.on_changed(update)

        plt.show()

        flagged_indices = np.array([edge1_slider.val, edge2_slider.val])

        template = pd.DataFrame(data={'energy_hist': energy_hist, 'bins': energy_bins[0:bin_plot_edge],
                                      'indices': np.append(flagged_indices,
                                                           np.zeros(len(energy_hist[0:bin_plot_edge]) - 2))})
        tl.make_path('Templates')
        template.to_csv(f'Templates/{self.unit}_{self.location["Location"]}_template.csv', index=False)
        if self.print_feedback:
            print('Template made')

    def _calibrate_NaI(self, energy_bins, energy_hist, spectra_conversions, spectra_frame):
        """Calibration algorithm for the sodium iodide scintillators."""
        flagged_indices = np.array([])
        # Takes the sum of each bin with its two closest neighboring bins on either side
        sums = energy_hist
        for i in range(2):
            sums += (np.roll(energy_hist, i + 1) + np.roll(energy_hist, i - 1))

        # Looks for the location of the maximum sum within the two bands where the peaks are likely to be
        band_starts = self.calibration_params['band_starts']
        band_ends = self.calibration_params['band_ends']
        for i in range(len(band_starts)):
            band_max = np.argmax(sums[band_starts[i]:band_ends[i]]) + int(band_starts[i])
            flagged_indices = np.append(flagged_indices, band_max)

        calibration = energy_bins[flagged_indices.astype(int)]
        if len(calibration) >= 2:
            print('For NaI:', file=spectra_conversions)
            print(f'{calibration[0]} V = 1.46 MeV', file=spectra_conversions)
            print(f'{calibration[1]} V = 2.60 MeV', file=spectra_conversions)

        spectra_frame['NaI'] = energy_hist
        self.update_attribute('NaI', 'calibration', calibration)
        return flagged_indices

    def _calibrate_LP(self, energy_bins, energy_hist, spectra_conversions, spectra_frame):
        """Calibration algorithm for the large plastic scintillators."""
        flagged_indices = np.array([])
        try:
            template = pd.read_csv(f'Templates/{self.unit}_{self.location["Location"]}_template.csv')
            correlation = signal.correlate(template['energy_hist'].to_numpy(), energy_hist, 'full')
            best_correlation_index = np.argmax(correlation)
            shift_amount = (-len(template) + 1) + best_correlation_index

            edge_indices = template['indices'].to_numpy()[0:2]
            flagged_indices = edge_indices + shift_amount
            self.lp_calibrated = True
        except FileNotFoundError:
            print('No LP template found for this location...', file=self.log)
            if self.print_feedback:
                print('No LP template found for this location...')

        calibration = energy_bins[flagged_indices.astype(int)]
        if len(calibration) >= 2:
            print('For LP:', file=spectra_conversions)
            print(f'{calibration[0]} V = 1.242 MeV', file=spectra_conversions)
            print(f'{calibration[1]} V = 2.381 MeV', file=spectra_conversions)

        spectra_frame['LP'] = energy_hist
        self.update_attribute('LP', 'calibration', calibration)
        return flagged_indices

    def _plot_spectra(self, scintillator, energy_bins, energy_hist, flagged_indices, sp_path):
        """Plots the histograms (with calibration lines, if applicable) for the given spectra."""
        # Plots the actual spectrum
        bin_plot_edge = len(energy_bins) - 1  # Histogram array is shorter than bin array by 1 (no idea why)
        bin_size = self.calibration_params['bin_size']
        plt.figure(figsize=[20, 11.0])
        plt.title(f'Energy Spectrum for {scintillator}, {self.full_date_str}', loc='center')
        plt.xlabel('Energy Channel')
        plt.ylabel('Counts/bin')
        plt.yscale('log')
        plt.bar(energy_bins[0:bin_plot_edge], energy_hist[0:bin_plot_edge], color='r',
                width=bin_size / 2, zorder=1)

        # Plots the energy bins corresponding to the desired energies as vertical lines
        if flagged_indices.size > 0:
            plt.vlines(energy_bins[flagged_indices.astype(int)], 0, np.amax(energy_hist), zorder=2, alpha=0.75)

        # Saves the figure
        plt.savefig(f'{sp_path}{scintillator}_Spectrum.png', dpi=500)
        plt.clf()

    def calibrate(self, existing_spectra=None):
        """Makes the energy spectra histograms and calibrates the large plastic and sodium iodide scintillators.
        Calibration energies for each scintillator are saved to their corresponding scintillator objects.

        Parameters
        ------
        existing_spectra : dict
            Optional. A dictionary whose entries correspond to energy spectra histograms for each scintillator.

        """

        # Fetching a few calibration parameters
        bin_range = self.calibration_params['bin_range']
        bin_size = self.calibration_params['bin_size']

        # Making the energy bins and setting up the calibration files
        energy_bins = np.arange(0.0, bin_range, bin_size)
        sp_path = f'{self.results_loc}Results/{self.unit}/{self.date_str}/'
        tl.make_path(sp_path)
        spectra_conversions = open(f'{sp_path}spectra_conversions.txt', 'w')
        spectra_frame = pd.DataFrame()
        spectra_frame['energy bins'] = energy_bins[:-1]
        if self or existing_spectra:
            if self.is_data_present('LP') or (existing_spectra and len(existing_spectra['LP']) != 0):
                energy_hist = self._generate_hist(energy_bins, 'LP', existing_spectra)  # Putting this up here
                # so that we don't have to do it again just for template mode
                if self.template:
                    self._make_template(energy_bins, energy_hist)

                # Calibrates the LP scintillator and plots the calibration
                flagged_indices = self._calibrate_LP(energy_bins, energy_hist, spectra_conversions, spectra_frame)
                self._plot_spectra('LP', energy_bins, energy_hist, flagged_indices, sp_path)
            else:
                print('Cannot calibrate LP (missing data)...', file=self.log)
                if self.print_feedback:
                    print('Cannot calibrate LP (missing data)...')

            # Calibrates the NaI scintillator (if possible) and plots the calibration
            if self.is_data_present('NaI') or (existing_spectra and len(existing_spectra['NaI']) != 0):
                energy_hist = self._generate_hist(energy_bins, 'NaI', existing_spectra)
                flagged_indices = self._calibrate_NaI(energy_bins, energy_hist, spectra_conversions, spectra_frame)
                self._plot_spectra('NaI', energy_bins, energy_hist, flagged_indices, sp_path)
            else:
                print('Cannot calibrate NaI (missing data)...', file=self.log)
                if self.print_feedback:
                    print('Cannot calibrate NaI (missing data)...')

            spectra_frame.to_json(f'{sp_path}{self.date_str}_spectra.json')
            spectra_conversions.close()

        else:
            raise ValueError("ValueError: data necessary for calibration is either missing or hasn't been imported.")