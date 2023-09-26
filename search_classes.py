"""Classes for use in the TGF and glow search program.

Classes:
    Detector:
        Object used to store all relevant information about the detector and the data.
    Chunk:
        Object used to store all relevant information about the detector and the data for a chunk of the day.
    ShortEvent:
        Object used to store all relevant information about potential short events.
    PotentialGlow:
        Object used to store all relevant information about potential glows.

"""

import glob as glob
import os as os
import contextlib as contextlib
import psutil as psutil
import numpy as np
import pandas as pd
import scipy.signal as signal
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import datetime as dt
import DataReaderFinal as dr
# import DataReaderTimetrack2 as dr
import search_module as sm


class Detector:
    """Used to store all relevant information about the detector and the data.

    The detector object is used to store the name of the detector, the requested date for analysis in various formats,
    and the actual data in a single, centralized location. Once data is imported from each file, it, along with the
    list of files, eRC serial number, and some other information, is stored in a nested dictionary structure that makes
    it very easy to access.

    Parameters
    ----------
    unit : str
        The name of the detector that the analysis is being requested for.
    first_sec : float
        The first second in EPOCH time of the day that data analysis is being requested for.
    modes : list
        A list of the requested modes that the program should operate under.

    Attributes
    ----------
    log : file
        The .txt file where program actions and findings are logged.
    full_day_string : str
        The timestamp for the requested day in yymmdd format.
    date_timestamp : str
        The timestamp for the requested in day in yyyy-mm-dd format.
    location : str
        The location of the detector on the requested day.
    import_path : str
        The directory path for the location of the requested data files.
    regex : function
        A lambda function that, when given the eRC serial number, returns the regex for a scintillator's files.
    good_lp_calibration : bool
        A flag for whether the program was able to calibrate the detector's large plastic scintillator or not.
    long_event_scint_list : list
        A list of the scintillators used by the long event search algorithm.
    calibration_params : dict
        A dictionary containing various parameters used in the detector calibration algorithm.
    scintillators : dict
        A nested dictionary containing all the relevant information for a detector's scintillators.
    THOR : bool
        A flag for whether the requested detector is a THOR unit or not.
    GODOT : bool
        A flag for whether the requested detector is GODOT or not.
    SANTIS : bool
        A flag for whether the requested detector is the Santis instrument or not.
    CROATIA : bool
        A flag for whether the requested detector is the Croatia instrument or not.
    custom : bool
        A flag for whether the program should operate in "custom" mode or not. Custom mode essentially just instructs
        the program to use a custom file path for importing data instead of the ones that are built-in. This path can be
        changed in the search_module module under the "C_raw_data" function.
    processed : bool
        A flag for whether the program should operate in "processed" mode or not. Under this mode, the program will use
        the built-in file path for GODOT processed data on Sol.
    template : bool
        A flag for whether the program should operate in "template" mode or not. This mode allows the user to generate
        and tweak templates that the program uses to calibrate the large plastic scintillator. These must be made for
        each new location.
    results_loc : str
        The directory where program results will be exported.
    GUI : bool
        A flag for whether the program is being executed from the GUI.

    """

    def __init__(self, unit, first_sec, modes):
        # Basic information
        self.unit = unit
        self.first_sec = first_sec
        self.log = None
        self.modes = modes
        self.full_day_string = dt.datetime.utcfromtimestamp(int(first_sec)).strftime('%y%m%d')  # In format yymmdd
        self.date_timestamp = dt.datetime.utcfromtimestamp(int(first_sec)).strftime('%Y-%m-%d')  # In format yyyy-mm-dd
        self.location = 'location'  # Eventually this will be fetched from a function in search_module
        self.good_lp_calibration = False

        # Detector information
        self.THOR = False
        self.GODOT = False
        self.SANTIS = False
        self.CROATIA = False

        if self.unit == 'GODOT':
            self.GODOT = True
            self.long_event_scint_list = ['NaI', 'LP']
            self.calibration_params = {'bin_range': 15008.0, 'bin_size': 16, 'band_starts': [38, 94],
                                       'band_ends': [75, 125], 'template_bin_plot_edge': 200}
            if 'processed' in self.modes:
                self.import_path = f'{sm.G_processed_data_loc()}/{self.full_day_string[0:4]}'
            else:
                self.import_path = f'{sm.G_raw_data_loc()}/{self.full_day_string}'

            self.regex = lambda eRC: f'eRC{eRC}_lm*_{self.full_day_string}_*'
            self.scintillators = {'NaI': {'eRC': '1490',
                                          'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                          'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                          'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                       'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7,
                                                       'started': 0}
                                          },
                                  'LP': {'eRC': '1491',
                                         'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                         'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                         'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                      'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7, 'started': 0}
                                         }
                                  }

        elif self.unit[0:4] == 'THOR':
            self.THOR = True
            self.long_event_scint_list = ['NaI']
            self.calibration_params = {'bin_range': 65535.0, 'bin_size': 1, 'band_starts': [2000, 5600],
                                       'band_ends': [2500, 6300], 'template_bin_plot_edge': 8000}  # placeholder
            self.import_path = f'{sm.T_raw_data_loc()}/{unit}/Data/{self.full_day_string}'
            self.regex = lambda eRC: f'eRC{eRC}*_lm_{self.full_day_string}_*'
            self.scintillators = {'NaI': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[0],
                                          'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                          'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                          'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                       'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7,
                                                       'started': 0}
                                          },
                                  'SP': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[1],
                                         'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                         'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                         'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                      'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7, 'started': 0}
                                         },
                                  'MP': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[2],
                                         'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                         'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                         'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                      'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7, 'started': 0}
                                         },
                                  'LP': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[3],
                                         'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                         'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                         'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                      'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7, 'started': 0}
                                         }
                                  }

        elif self.unit == 'SANTIS':
            self.SANTIS = True
            self.long_event_scint_list = ['LP']
            self.calibration_params = {'bin_range': 15008.0, 'bin_size': 16, 'band_starts': [38, 94],
                                       'band_ends': [75, 125], 'template_bin_plot_edge': 400}
            self.import_path = f'{sm.S_raw_data_loc()}/{self.full_day_string}'
            self.regex = lambda eRC: f'eRC{eRC}*_lm_{self.full_day_string}_*'
            self.scintillators = {'LP': {'eRC': '2549',
                                         'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                         'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                         'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                      'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7, 'started': 0}
                                         }
                                  }

        elif self.unit == 'CROATIA':
            self.CROATIA = True
            self.long_event_scint_list = ['LP']
            self.calibration_params = {'bin_range': 15008.0, 'bin_size': 16, 'band_starts': [38, 94],
                                       'band_ends': [75, 125], 'template_bin_plot_edge': 400}
            self.import_path = f'{sm.CR_raw_data_loc()}/{self.full_day_string}'
            self.regex = lambda eRC: f'eRC{eRC}*_lm_{self.full_day_string}_*'
            self.scintillators = {'MP': {'eRC': '4193',
                                         'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                         'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                         'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                      'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7, 'started': 0}
                                         },
                                  'LP': {'eRC': '2549',
                                         'filelist': [], 'filetime_extrema': [], 'calibration': [],
                                         'time': np.array([]), 'energy': np.array([]), 'wc': np.array([]),
                                         'passtime': {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0,
                                                      'ppsunix': -1.0, 'lastwc': 0, 'ppswc': 0, 'hz': 8e7, 'started': 0}
                                         },
                                  }

        else:
            raise ValueError('ValueError: not a valid detector.')

        # Modes
        self.custom = False
        self.processed = False

        if 'custom' in self.modes:
            self.custom = True  # Not necessary for anything right now
            self.import_path = f'{sm.C_raw_data_loc()}'

        if 'processed' in self.modes:
            self.processed = True
            if self.unit != 'GODOT':
                raise ValueError('ValueError: processed data is only accessible for GODOT')

        self.template = True if 'template' in self.modes else False

        self.results_loc = sm.results_loc()
        self.GUI = True if 'GUI' in self.modes else False
        if self.GUI:
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

    # String casting magic method
    def __str__(self):
        return f'Detector({self.unit}, {self.first_sec}, {self.modes})'

    # Debugging string magic method
    def __repr__(self):
        data_attributes = ['time', 'energy', 'wc']
        scintillators_with_data = []
        has_data = False
        for scintillator in self.scintillators:
            medium = self.scintillators[scintillator]
            for attribute in data_attributes:
                if len(medium[attribute]) > 0:
                    has_data = True
                    scintillators_with_data.append(scintillator)
                    break

        default_string = self.__str__()
        data_string = f' in {scintillators_with_data}' if has_data else ''
        return default_string + f' Has data = {has_data}' + data_string

    # Iterator for use in loops
    def __iter__(self):
        scintillator_keys = list(self.scintillators.keys())
        for scintillator in scintillator_keys:
            yield scintillator

    # Bool casting overload. Returns True if data for necessary scintillators (LP, NaI) is present
    def __bool__(self):
        if 'LP' in self.long_event_scint_list:
            necessary_scintillators = self.long_event_scint_list
        else:
            necessary_scintillators = self.long_event_scint_list + ['LP']

        data_present = True
        for scint in necessary_scintillators:
            # Using energy as the check here is arbitrary; It could just as easily have been time instead
            if len(self.attribute_retriever(scint, 'energy')) == 0:
                data_present = False

        return data_present

    def attribute_retriever(self, scintillator, attribute):
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
            if attribute in medium:
                desired_attribute = medium[attribute]
                return desired_attribute
            else:
                raise AttributeError('AttributeError: not a valid attribute')

        else:
            raise AttributeError('AttributeError: not a valid scintillator')

    def attribute_updator(self, scintillator, attribute, new_info):
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
                    raise Exception('Attribute and new info must be the same length')

                for i in range(len(attribute)):
                    if attribute[i] in medium:
                        medium.update({attribute[i]: new_info[i]})
                    else:
                        raise AttributeError('AttributeError: not a valid attribute')

            elif type(attribute) is str:
                if attribute in medium:
                    medium.update({attribute: new_info})
                else:
                    raise AttributeError('AttributeError: not a valid attribute')

            else:
                raise TypeError('TypeError: if attribute is a list, new_info must also be a list')

            self.scintillators.update({scintillator: medium})

        else:
            raise AttributeError('AttributeError: not a valid scintillator')

    def _generate_hist(self, energy_bins, scintillator, existing_spectra=None):
        """Returns the energy spectra histogram for the requested scintillator."""
        energies = self.attribute_retriever(scintillator, 'energy')
        if existing_spectra is not None:
            energy_hist = existing_spectra[scintillator]
        else:
            energy_hist, bin_edges = np.histogram(energies, bins=energy_bins)

        return energy_hist

    def _make_template(self, energy_bins, energy_hist):
        """Makes a template that can be used in the LP calibration algorithm's cross-correlation."""
        bin_plot_edge = len(energy_bins) - 1  # Histogram array is shorter than bin array by 1 (no idea why)

        template_bin_plot_edge = self.calibration_params['template_bin_plot_edge']

        sm.print_logger('Entering template mode...', self.log)
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
        sm.path_maker('Templates')
        template.to_csv(f'Templates/{self.unit}_{self.location}_template.csv', index=False)
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
        self.attribute_updator('NaI', 'calibration', calibration)
        return flagged_indices

    def _calibrate_LP(self, energy_bins, energy_hist, spectra_conversions, spectra_frame):
        """Calibration algorithm for the large plastic scintillators."""
        flagged_indices = np.array([])
        try:
            template = pd.read_csv(f'Templates/{self.unit}_{self.location}_template.csv')
            correlation = signal.correlate(template['energy_hist'].to_numpy(), energy_hist, 'full')
            best_correlation_index = np.argmax(correlation)
            shift_amount = (-len(template) + 1) + best_correlation_index

            edge_indices = template['indices'].to_numpy()[0:2]
            flagged_indices = edge_indices + shift_amount
            self.good_lp_calibration = True
        except FileNotFoundError:
            sm.print_logger('No LP template found for this location...', self.log)

        calibration = energy_bins[flagged_indices.astype(int)]
        if len(calibration) >= 2:
            print('For LP:', file=spectra_conversions)
            print(f'{calibration[0]} V = 1.242 MeV', file=spectra_conversions)
            print(f'{calibration[1]} V = 2.381 MeV', file=spectra_conversions)

        spectra_frame['LP'] = energy_hist
        self.attribute_updator('LP', 'calibration', calibration)
        return flagged_indices

    def _plot_spectra(self, scintillator, energy_bins, energy_hist, flagged_indices, sp_path):
        """Plots the histograms (with calibration lines, if applicable) for the given spectra."""
        # Plots the actual spectrum
        bin_plot_edge = len(energy_bins) - 1  # Histogram array is shorter than bin array by 1 (no idea why)
        bin_size = self.calibration_params['bin_size']
        plt.figure(figsize=[20, 11.0])
        plt.title(f'Energy Spectrum for {scintillator}, {self.date_timestamp}', loc='center')
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
        Calibration energies for each scintillator are saved to the 'calibration' attribute of the detector object's
        nested dictionary structure.

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
        sp_path = f'{sm.results_loc()}Results/{self.unit}/{self.full_day_string}/'
        sm.path_maker(sp_path)
        spectra_conversions = open(f'{sp_path}spectra_conversions.txt', 'w')
        spectra_frame = pd.DataFrame()
        spectra_frame['energy bins'] = energy_bins[:-1]
        if self or existing_spectra:
            energy_hist = self._generate_hist(energy_bins, 'LP', existing_spectra)  # Putting this up here
            # so that we don't have to do it again just for template mode
            if self.template:
                self._make_template(energy_bins, energy_hist)

            # Calibrates the LP scintillator and plots the calibration
            flagged_indices = self._calibrate_LP(energy_bins, energy_hist, spectra_conversions, spectra_frame)
            self._plot_spectra('LP', energy_bins, energy_hist, flagged_indices, sp_path)

            # Calibrates the NaI scintillator and plots the calibration
            energy_hist = self._generate_hist(energy_bins, 'NaI', existing_spectra)
            flagged_indices = self._calibrate_NaI(energy_bins, energy_hist, spectra_conversions, spectra_frame)
            self._plot_spectra('NaI', energy_bins, energy_hist, flagged_indices, sp_path)

            spectra_frame.to_json(f'{sp_path}{self.full_day_string}_spectra.json')
            spectra_conversions.close()
        else:
            raise ValueError("ValueError: Data for calibration is either missing or hasn't been imported")

    def _filter_files(self, complete_filelist):
        """Returns an ordered list of files with duplicate/incompatible files filtered out."""
        # Filters out trace mode files and .txtp files (whatever those are)
        filtered_filelist = []
        unfiltered_filelist = []
        for j in range(len(complete_filelist)):
            current_file = complete_filelist[j]
            if current_file[-3:] == 'xtr' or current_file[-4:] == 'txtp' or current_file[-5:] == 'xtrpp':
                continue
            elif current_file[-7:] == '.txt.gz':
                filtered_filelist.append(current_file)
            else:
                unfiltered_filelist.append(current_file)

        # Eliminates duplicate files
        for j in range(len(unfiltered_filelist)):
            current_file = unfiltered_filelist[j]
            if f'{current_file}.gz' in filtered_filelist:
                continue
            else:
                filtered_filelist.append(current_file)

        # Puts the files in the correct order
        files = []
        extensions = []
        for j in range(len(filtered_filelist)):
            file = filtered_filelist[j]
            if file[-4:] == '.txt':
                files.append(file.replace('.txt', ''))
                extensions.append('.txt')
            elif file[-4:] == '.csv':
                files.append(file.replace('.csv', ''))
                extensions.append('.csv')
            elif file[-7:] == '.txt.gz':
                files.append(file.replace('.txt.gz', ''))
                extensions.append('.txt.gz')

        file_order = np.argsort(files)
        files.sort()
        extensions = [extensions[s] for s in file_order]
        filelist = []
        for j in range(len(files)):
            filelist.append(f'{files[j]}{extensions[j]}')

        return filelist

    def calculate_fileset_size(self):
        """Returns the total size (in bytes) of all the currently held files for the day."""
        total_file_size = 0
        for scintillator in self:
            filelist = self.attribute_retriever(scintillator, 'filelist')
            for file in filelist:
                total_file_size += os.path.getsize(file)

        return total_file_size

    def data_importer(self, existing_filelists=False):
        """Imports data from data files into arrays and then updates them into the detector object's
        nested dictionary structure.

        Parameters
        ----------
        existing_filelists : bool
            Optional. If True, the function will use the file lists already stored in the detector object's nested
            dictionary structure.

        """

        for scintillator in self:
            if existing_filelists:
                break

            eRC = self.attribute_retriever(scintillator, 'eRC')
            # Here in case the data files in a custom location are grouped into daily folders
            try:
                complete_filelist = glob.glob(f'{self.import_path}/{self.regex(eRC)}')
                assert len(complete_filelist) > 0, 'Empty filelist'

            except AssertionError:
                complete_filelist = glob.glob(f'{self.import_path}/{self.full_day_string}'
                                              f'/{self.regex(eRC)}')

            filelist = self._filter_files(complete_filelist)
            self.attribute_updator(scintillator, 'filelist', filelist)

        # Checks to see if the necessary files for a full search are present
        if 'LP' in self.long_event_scint_list:
            necessary_scintillators = self.long_event_scint_list
        else:
            necessary_scintillators = self.long_event_scint_list + ['LP']

        data_present = True
        missing_data_scints = []
        for scint in necessary_scintillators:
            if len(self.attribute_retriever(scint, 'filelist')) == 0:
                data_present = False
                missing_data_scints.append(scint)

        if not data_present:
            print('\n')
            print('\n', file=self.log)
            sm.print_logger('No/missing necessary data for specified day.', self.log)
            print(f'Necessary data missing in the following: {", ".join(missing_data_scints)}', file=self.log)
            print('\n')
            raise FileNotFoundError

        # Determines whether there is enough free memory to load the entire dataset
        total_file_size = self.calculate_fileset_size()
        operating_memory = sm.memory_allowance()
        available_memory = psutil.virtual_memory()[1]/4
        if (operating_memory + total_file_size) > available_memory:
            raise MemoryError('MemoryError: not enough free memory to hold complete dataset.')

        for scintillator in self.scintillators:
            eRC = self.attribute_retriever(scintillator, 'eRC')
            filelist = self.attribute_retriever(scintillator, 'filelist')
            sm.print_logger('\n', self.log)
            sm.print_logger(f'For eRC {eRC} ({scintillator}):', self.log)

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
                    if not self.GUI:
                        print(f'{files_imported}/{len(filelist)} files imported', end='\r')
                    elif self.GUI and filecount_switch:
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
                                # data, passtime = dr.fileNameToData(file,
                                #                                    self.attribute_retriever(scintillator, 'passtime'))
                                data = dr.fileNameToData(file)

                            # self.attribute_updator(scintillator, 'passtime', passtime)
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

                print(f'{files_imported}/{len(filelist)} files imported', end='\r')

                # Makes the final arrays and exports them
                times = np.concatenate(time_list)
                # Corrects for the fact that the first 200-300 seconds of the next day are included in the last file
                day_change_array = np.array(np.where(np.diff(times) < -80000))
                if day_change_array.size > 0:
                    change_index = int(day_change_array[0]) + 1
                    times = np.append(times[0:change_index], times[change_index:] + 86400.0)

                # Does it for the file time extrema too
                for k in range(int(len(filetime_extrema_list)/8)):  # Last eighth of the files
                    last_file_extrema = filetime_extrema_list[-(k+1)]
                    for j in range(2):
                        extrema = last_file_extrema[j]
                        if extrema < 500:  # Extrema belonging to the next day will always be < 500
                            last_file_extrema[j] = extrema + 86400

                    filetime_extrema_list[-(k+1)] = last_file_extrema

                updated_attributes = ['time', 'energy', 'wc', 'filetime_extrema']
                updated_info = [times, np.concatenate(energy_list), np.concatenate(wallclock_list),
                                filetime_extrema_list]
                self.attribute_updator(scintillator, updated_attributes, updated_info)

                print('\n', file=self.log)
                print(f'Total Counts: {len(np.concatenate(time_list))}', file=self.log)
                print(f'Average time gap: {np.sum(file_time_gaps) / len(file_time_gaps)}', file=self.log)
                print('\n', file=self.log)

            except AssertionError:
                print('Missing data for the specified day.', file=self.log)
                print('Missing data for the specified day.', end='\r')
                continue

            except Exception as ex:
                if str(ex) == 'Reader Error':
                    print('Error with data reader.', file=self.log)
                    print('Error with data reader.', end='\r')
                    continue
                else:
                    raise


class Chunk(Detector):
    """Used to store all relevant information about the detector and the data for a chunk of the day.

    This is a child class of Detector, see Detector documentation for full list of methods and attributes. Chunk is used
    to store information and data for a particular chunk of the day. It is meant to be used in search.py's low memory
    mode in place of regular Detector objects, and as such it contains a few extra features.


    Parameters
    ----------
    unit : str
        The name of the particular detector that the analysis is being requested for.
    first_sec : float
        The first second in EPOCH time of the day that data analysis is being requested for.
    modes : list
        A list of the requested modes that the program should operate under.

    Attributes
    ----------
    scint_list : list
        A list of the detector's scintillators.

    """

    def __init__(self, unit, first_sec, modes):
        super().__init__(unit, first_sec, modes)
        self.scint_list = list(self.scintillators.keys())

    def return_passtime(self):
        """Returns a dictionary containing the passtime attributes for each of the detector's scintillators.

        Returns
        ----------
        dict
            A dictionary containing the passtime dictionaries for each of the detector's scintillators.

        """

        passtime_dict = {}
        for scintillator in self:
            passtime = self.attribute_retriever(scintillator, 'passtime')
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
            self.attribute_updator(scintillator, 'passtime', passtime_dict[scintillator])

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
            energies = self.attribute_retriever(scintillator, 'energy')
            chunk_hist, bin_edges = np.histogram(energies, bins=energy_bins)
            if len(existing_spectra_dict[scintillator]) == 0:
                existing_spectra_dict.update({scintillator: chunk_hist})
            else:
                new_hist = existing_spectra_dict[scintillator] + chunk_hist
                existing_spectra_dict.update({scintillator: new_hist})

        return existing_spectra_dict


class ShortEvent:
    """Object used to store all relevant information about potential short events.

    Parameters
    ----------
    event_start : int
        The index of the time array which corresponds to the start of the event.
    event_length : int
        The number of entries in the time array which make up the event.
    scintillator : str
        A string corresponding to the scintillator which the event was found in.

    Attributes
    ----------
    start : int
        The index of the time array which corresponds to the start of the event.
    length : int
        The number of entries in the time array which make up the event.
    stop : int
        The index of the time array which corresponds to the end of the event.
    scintillator : str
        A string corresponding to the scintillator which the event was found in.

    """

    def __init__(self, event_start, event_length, scintillator):
        self.start = int(event_start)
        self.length = int(event_length)
        self.stop = int(event_start + event_length)
        self.scintillator = scintillator

    # String casting magic method
    def __str__(self):
        return f'{self.scintillator} short event at index {self.start}'

    # Debugging string magic method
    def __repr__(self):
        return f'{self.scintillator} short event; start:{self.start}; stop:{self.stop}; length:{self.length}'

    def get_filename(self, times, filelist, filetime_extrema):
        """Gets the name of the file that the event occurred in.

        Parameters
        ----------
        times : np.array
            A numpy array containing times for each count.
        filelist : list
            A list of all the files to be searched when looking for the beginning of an event. Once the name of the
            file for an event has been found it is added to the scatter plot title.
        filetime_extrema : list
            A list of arrays containing the first and last second in each file.

        Returns
        -------
        filename : str
            The name of the file that the event occurred in.
        new_filelist : list
            A shorter version of the filelist for the requested scintillator. This new list is
            eventually used when the next scatter plot is generated to make file finding faster.
        new_filetime_extrema : list
            A shorter version of the file extrema for the requested scintillator. This new list is
            eventually used when the next scatter plot is generated to make file finding faster.

        """

        event_time = times[self.start] - 86400 if times[self.start] > 86400 else times[self.start]

        event_file = ''
        files_examined = 0
        new_filelist = []
        new_filetime_extrema = []
        for i in range(len(filetime_extrema)):
            first_time = filetime_extrema[i][0]
            last_time = filetime_extrema[i][1]
            if first_time <= event_time < last_time:
                event_file = filelist[i]
                new_filelist = filelist[files_examined:]
                new_filetime_extrema = filetime_extrema[files_examined:]
                break

            files_examined += 1

        return event_file, new_filelist, new_filetime_extrema

    def json_maker(self, detector, times, energies, wallclock, event_number, event_file):
        """Makes the short event JSON files.

        Parameters
        ----------
        detector : Detector object
            The detector object used to store all the data and relevant information.
        times : np.array
            A numpy array containing times for each count.
        energies : np.array
            A numpy array containing energies for each count.
        wallclock : np.array
            A numpy array containing wallclock times for each count.
        event_number : int
            A number corresponding to the event's number for that particular scintillator (i.e. 1 would be the first
            event for whatever scintillator, 2 would be the second, and so on).
        event_file : str
            The name of the file that the event occurred in.

        """

        event_times = times[self.start:self.stop]
        event_energies = energies[self.start:self.stop]
        event_wallclock = wallclock[self.start:self.stop]

        eventpath = (f'{detector.results_loc}Results/{detector.unit}/'
                     f'{detector.full_day_string}/event files/short events/')
        sm.path_maker(eventpath)
        event_frame = pd.DataFrame()
        event_frame['wc'] = event_wallclock
        event_frame['SecondsOfDay'] = event_times
        event_frame['energies'] = event_energies
        event_frame['file'] = event_file  # Note: this column will be filled by the same file name over and over again

        # Saves the json file
        event_frame.to_json(f'{eventpath}{detector.full_day_string}_{self.scintillator}_event{event_number}.json')

    def scatterplot_maker(self, timescales, detector, times, energies, event_number, event_file):
        """Makes the short event scatter plots.

        Parameters
        ----------
        timescales : list
            A list of the timescales (in seconds) that the scatter plots are generated in.
        detector : Detector object
            The detector object used to store all the data and relevant information.
        times : np.array
            A numpy array containing times for each count.
        energies : np.array
            A numpy array containing energies for each count.
        event_number : int
            A number corresponding to the event's number for that particular scintillator (i.e. 1 would be the first
            event for whatever scintillator, 2 would be the second, and so on).
        event_file : str
            The name of the file that the event occurred in.

        """

        # Truncated time and energy arrays to speed up scatter plot making
        fraction_of_day = 1/64
        spacer = int((len(times)*fraction_of_day)/2)
        left_edge = 0 if self.start-spacer < 0 else self.start - spacer
        right_edge = (len(times) - 1) if self.stop + spacer > (len(times) - 1) else self.stop + spacer
        trunc_times = times[left_edge:right_edge]
        trunc_energies = energies[left_edge:right_edge]

        event_times = times[self.start:self.stop]
        event_energies = energies[self.start:self.stop]
        event_length = event_times[-1] - event_times[0]

        figure1 = plt.figure(figsize=[20, 11.0])
        figure1.suptitle(f'{self.scintillator} Event {str(event_number)}, '
                         f'{dt.datetime.utcfromtimestamp(times[self.start] + detector.first_sec)} UTC, '
                         f'{len(event_energies)} counts', fontsize=20)
        ax1 = figure1.add_subplot(3, 1, 1)
        ax2 = figure1.add_subplot(3, 1, 2)
        ax3 = figure1.add_subplot(3, 1, 3)
        ax_list = [ax1, ax2, ax3]

        for i in range(len(ax_list)):
            ts = timescales[i]
            ax = ax_list[i]
            padding = (ts - event_length) / 2
            if event_length >= ts:
                best_time = event_times[np.argmin(np.abs(event_times - np.roll(event_times, 1)))]
                ax.set_xlim(xmin=best_time - (ts / 2), xmax=best_time + (ts / 2))
            else:
                ax.set_xlim(xmin=event_times[0] - padding, xmax=event_times[-1] + padding)

            dot_size = 3 if ts == timescales[0] else 1  # makes larger dots for top plot
            ax.set_yscale('log')
            ax.set_ylim([0.5, 1e5])
            ax.scatter(trunc_times, trunc_energies + 0.6, s=dot_size, zorder=1, alpha=1.0)
            ax.set_xlabel(f'Time (Seconds, {ts}s total)')
            ax.set_ylabel('Energy Channel')
            # Lines appear (100*percent)% to the left or right of event start/stop depending on subplot timescale
            percent = 0.001
            # Event start line is orange, end line is red
            ax.vlines([event_times[0] - percent*ts, event_times[-1] + percent*ts], 0, 1e5,
                      colors=['orange', 'r'], linewidth=1, zorder=-1, alpha=0.3)

        # Adds the name of the relevant data file to the scatter plot
        plt.title(f'Obtained from {event_file}', fontsize=15, y=-0.4)

        # Saves the scatter plot
        # Note: with this code, if an event happens in that 200-300 seconds of the next day that are included in the
        # last file, the image will have the wrong date in its name (though the timestamp in the scatter plot title will
        # always be correct)
        scatterpath = (f'{detector.results_loc}Results/{detector.unit}/'
                       f'{detector.full_day_string}/scatterplots/')
        sm.path_maker(scatterpath)
        figure1.savefig(f'{scatterpath}{detector.full_day_string}_{self.scintillator}_event{event_number}.png')
        plt.close(figure1)


class PotentialGlow:
    """Object used to store all relevant information about potential glows.

    Parameters
    ----------
    glow_start : int
        The index of the histogram bin which corresponds to the start of the event.
    glow_length : int
        The number of histogram bins which make up an event.

    Attributes
    ----------
    start : int
        The index of the histogram bin which corresponds to the start of the event.
    length : int
        The number of histogram bins which make up an event.
    stop : int
        The index of the histogram bin which corresponds to the end of the event.
    peak_index : int
        The location of the bin with the largest z-score among all the bins comprising the event.
    highest_score : float
        The largest z-score in the event.
    start_sec : int
        The beginning of the event in seconds.
    stop_sec : int
        The end of the event in seconds.

    """

    def __init__(self, glow_start, glow_length):
        self.start = int(glow_start)
        self.length = int(glow_length)
        self.stop = int(glow_start + glow_length - 1) if self.length > 1 else int(glow_start + glow_length)
        self.peak_index = 0
        self.highest_score = 0
        self.start_sec = 0
        self.stop_sec = 0

    # String casting magic method
    def __str__(self):
        return f'Long event at index {self.start}'

    # Debugging string magic method
    def __repr__(self):
        return f'Long event; start:{self.start}; stop:{self.stop}; length:{self.length}'

    def highest_zscore(self, z_scores):
        """Identifies the highest z-score and its corresponding bin for an event."""
        glow_scores = z_scores[self.start:self.stop]
        highest_score = np.max(glow_scores)
        self.peak_index = np.argmax(glow_scores) + self.start
        return highest_score

    def beginning_and_end_seconds(self, day_bins, binsize):
        """Retrieves the beginning and total length of an event in seconds."""
        glow_times = day_bins[self.start:self.stop]
        first_sec = glow_times[0]
        length = self.length * binsize
        last_sec = first_sec + length
        return first_sec, last_sec

    def hist_subplotter(self, ax, day_bins, hist_allday, mue, sigma, flag_threshold):
        """Makes the flagged z-score subplots for the glow search algorithm's full day histogram.

        The subplots generated by this function are meant to highlight the four most interesting-looking potential glows
        for a given day.

        Parameters
        ----------
        ax : plt.ax
            The pyplot axis object that the subplot is going to be generated for.
        day_bins : np.array
            The numpy array specifying the bins used to make the histogram.
        hist_allday : np.array
            The numpy array containing the counts per time bin.
        mue : np.array
            The average count rates expected at each bin (constant for non-aircraft mode).
        sigma : np.array
            The standard deviation of each bin's expected average (constant for non-aircraft mode).
        flag_threshold : int
            The z-score above mue threshold at which glows are flagged.

        """
        padding = 20
        c = 0 if (self.peak_index - padding) < 0 else (self.peak_index - padding)
        d = (len(day_bins) - 2) if (self.peak_index + padding) > (len(day_bins) - 2) else (self.peak_index + padding)

        subbins = day_bins[c:d]
        subhist = hist_allday[c:d]
        submue = mue[c:d]
        subsigma = sigma[c:d]
        binsize = int(day_bins[1] - day_bins[0])
        ax.bar(subbins, subhist, alpha=0.5, color='c', width=binsize)
        ax.set_xlabel('Seconds of Day (UT)')
        ax.set_ylabel('Counts/bin')
        ax.plot(subbins, submue + flag_threshold * subsigma, color='blue', linestyle='dashed', linewidth=2)
        ax.grid(True)
