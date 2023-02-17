"""Classes for use in the TGF and glow search program.

Classes:
    Detector:
        Object used to store all relevant information about the detector and the data.
    ShortEvent:
        Object used to store all relevant information about potential short events.
    PotentialGlow:
        Object used to store all relevant information about potential glows.

"""

import glob
import numpy as np
import pandas as pd
import scipy.signal as signal
import matplotlib.pyplot as plt
import search_module as sm
import DataReaderFinal as dr
import datetime as dt


class Detector:
    """ Used to store all relevant information about the detector and the data.

    The detector object is used to store the name of the detector, the requested date for analysis in various formats,
    and the actual data in a single, centralized location. Once data is imported from each file, it, along with the
    list of files and the eRC serial number, is stored in a nested dictionary structure that makes it very easy to
    access.

    Parameters
    ----------
    unit : str
        The name of the particular detector that the analysis is being requested for.
    first_sec : float
        The first second in EPOCH time of the day that data analysis is being requested for.
    log : file
        The .txt file where program actions and findings are logged.
    modes : list
        A list of the requested modes that the program should operate under.

    Attributes
    ----------
    full_day_string : str
        The timestamp for the requested day in yymmdd format.
    date_timestamp : str
        The timestamp for the requested in day in yyyy-mm-dd format.
    location : str
        The location of the detector on the requested day.
    import_path : str
        The directory path for the location of the requested data files.
    good_lp_calibration : bool
        A flag for whether the program was able to calibrate the detector's large plastic scintillator or not.
    THOR : bool
        A flag for whether the requested detector is a THOR unit or not.
    GODOT : bool
        A flag for whether the requested detector is GODOT or not.
    SANTIS : bool
        A flag for whether the requested detector is the Santis instrument or not.
    custom : bool
        A flag for whether the program should operate in "custom" mode or not. Custom mode essentially just instructs
        the program to use a custom file path for importing data instead of the ones that are built-in. This path can be
        changed in the search_module module under the "C_raw_data" function.
    processed : bool
        A flag for whether the program should operate in "processed" mode or not. Under this mode, the program will use
        the built-in file path for GODOT processed data on Sol.
    plastics : bool
        A flag for whether the program should operate in "plastics" mode or not. When this mode is requested, the
        program will run the short event search algorithm on the small, medium, and large plastic scintillator data
        as opposed to just the large plastic.
    template : bool
        A flag for whether the program should operate in "template" mode or not. This mode allows the user to generate
        and tweak templates that the program uses to calibrate the large plastic scintillator. These must be made for
        each new location.

    """
    def __init__(self, unit, first_sec, log, modes):
        # Basic information
        self.unit = unit
        self.first_sec = first_sec
        self.log = log
        self.modes = modes
        self.full_day_string = dt.datetime.utcfromtimestamp(int(first_sec)).strftime('%y%m%d')  # In format yymmdd
        self.date_timestamp = dt.datetime.utcfromtimestamp(int(first_sec)).strftime('%Y-%m-%d')  # In format yyyy-mm-dd
        self.location = 'location'  # Eventually this will be fetched from a function in search_module
        self.import_path = ''
        self.good_lp_calibration = False

        # Detector information
        self.THOR = False
        self.GODOT = False
        self.SANTIS = False

        if self.unit == 'GODOT':
            self.GODOT = True
            if 'processed' in self.modes:
                self.import_path = f'{sm.G_processed_data_loc()}/{self.full_day_string[0:4]}'
            else:
                self.import_path = f'{sm.G_raw_data_loc()}/{self.full_day_string}'

            self.scintillators = {'NaI': {'eRC': '1490', 'filelist': [],
                                          'time': np.array([]), 'energy': np.array([])},
                                  'LP': {'eRC': '1491', 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}

        elif self.unit[0:4] == 'THOR':
            self.THOR = True
            self.import_path = f'{sm.T_raw_data_loc()}/{unit}/Data/{self.full_day_string}'
            self.scintillators = {'NaI': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[0], 'filelist': [],
                                          'time': np.array([]), 'energy': np.array([])},
                                  'SP': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[1], 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])},
                                  'MP': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[2], 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])},
                                  'LP': {'eRC': sm.T_eRC(self.unit, self.full_day_string)[3], 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}

        elif self.unit == 'SANTIS':
            self.SANTIS = True
            self.import_path = f'{sm.S_raw_data_loc()}/{self.full_day_string}'
            self.scintillators = {'LP': {'eRC': '2549', 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}
        else:
            print('Not a valid detector')
            exit()

        # Modes
        self.custom = False
        self.processed = False
        self.plastics = False
        self.template = False

        if 'custom' in self.modes:
            self.custom = True  # Not necessary for anything right now
            self.import_path = f'{sm.C_raw_data_loc()}'

        if 'processed' in self.modes:
            self.processed = True
            if self.unit != 'GODOT':
                print('Processed data is only accessible for GODOT')
                exit()

        if 'plastics' in self.modes:
            self.plastics = True

        if 'template' in self.modes:
            self.template = True

    def attribute_retriever(self, scintillator, attribute):  # Make it possible to request multiple attributes at once?
        """Retrieves the requested attribute for a particular scintillator.

        Parameters
        ----------
        scintillator : str
            A string corresponding to the scintillator of interest. Allowed values: 'NaI', 'SP', 'MP', 'LP'.
        attribute : str
            A string corresponding to the scintillator attribute of interest. Allowed values: 'eRC', 'filelist', 'time',
            'energy'.

        Returns
        -------
        str, list, np.array
            String if 'eRC' is requested, list if 'filelist' is requested, or a numpy array full of float values if
            'time' or 'energy' is requested.

        """

        medium = self.scintillators[scintillator]
        desired_attribute = medium[attribute]
        return desired_attribute

    def spectra_maker(self):
        """Makes the energy spectra histograms and calibrates the large plastic and sodium iodide scintillators.

        Returns
        -------
        np.array
            Two numpy arrays. The first array contains the large plastic scintillator energy channels corresponding to
            the most compton-scattered photons of Potassium-40 and Thorium (i.e. compton edge locations). The second
            array contains the sodium iodide scintillator energy channels corresponding to the photo peaks of
            Potassium-40 and Thorium.

        """
        lp_energies = np.array([])
        nai_energies = np.array([])
        if self.THOR:
            bin_range = 65535.0
            bin_number = 65536
            band_starts = [2000, 5600]
            band_ends = [2500, 6300]
            template_bin_plot_edge = 8000
        else:
            bin_range = 15008.0
            bin_number = 939
            band_starts = [38, 94]
            band_ends = [75, 125]
            template_bin_plot_edge = 200

        for scintillator in self.scintillators:
            if scintillator == 'SP' or scintillator == 'MP':
                continue  # (for now)

            energies = self.attribute_retriever(scintillator, 'energy')
            energy_bins = np.linspace(0.0, bin_range, num=bin_number)
            energy_hist, bin_edges = np.histogram(energies, bins=energy_bins)
            bin_plot_edge = len(energy_bins)
            flagged_indices = np.array([])
            # Makes a template that can be used in the LP calibration algorithm's cross-correlation
            if self.template and scintillator == 'LP':
                sm.print_logger('Entering template mode...', self.log)
                print('\n')
                iterations = 0
                bin_plot_edge = template_bin_plot_edge
                while True:
                    if iterations > 0:
                        print(f'Previous line locations: K40: {int(flagged_indices[0])}, T: {int(flagged_indices[1])}')
                    else:
                        if self.THOR:
                            print('0 and 0 are good starting positions for THOR')  # fix this
                        else:
                            print('41 and 83 are good starting positions for GODOT')

                    flagged_indices = np.array([])
                    edge1 = 0
                    edge2 = 0
                    while True:
                        try:  # Idiot proofing
                            edge1 = int(input('K40 Line Location: '))
                            edge2 = int(input('T Line Location: '))
                            break
                        except ValueError:
                            print('Only one NUMBER at a time please')
                            continue

                    flagged_indices = np.append(flagged_indices, edge1)
                    flagged_indices = np.append(flagged_indices, edge2)

                    plt.xlabel('Energy Channel')
                    plt.ylabel('Counts/bin')
                    plt.yscale('log')
                    plt.hist(energies, bins=energy_bins[0:bin_plot_edge], color='r', rwidth=0.5, zorder=1)
                    plt.vlines(energy_bins[flagged_indices.astype(int)], 0, 1e6, zorder=2, alpha=0.75)
                    plt.show()

                    adequate = input('Are these good locations? (y/n): ')
                    print('\n')
                    if adequate in ['Y', 'y']:
                        break
                    iterations += 1

                template = pd.DataFrame(data={'energy_hist': energy_hist, 'bins': energy_bins[0:938],
                                              'indices': np.append(flagged_indices, np.zeros(len(energy_hist)-2))})
                sm.path_maker('Templates')
                template.to_csv(f'Templates/{self.unit}_{self.location}_template.csv', index=False)
                print('Template made')
            elif scintillator == 'LP':
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

                lp_energies = energy_bins[flagged_indices.astype(int)]

            else:  # NaI scintillator
                # Takes the sum of each bin with its two closest neighboring bins on either side
                sums = energy_hist
                for i in range(2):
                    sums += (np.roll(energy_hist, i + 1) + np.roll(energy_hist, i - 1))

                # Looks for the location of the maximum sum within the two bands where the peaks are likely to be
                for i in range(len(band_starts)):
                    band_max = np.argmax(sums[band_starts[i]:band_ends[i]]) + int(band_starts[i])
                    flagged_indices = np.append(flagged_indices, band_max)
                nai_energies = energy_bins[flagged_indices.astype(int)]

            # Plots the actual spectrum
            plt.figure(figsize=[20, 11.0])
            plt.title(f'Energy Spectrum for {scintillator}, {self.date_timestamp}', loc='center')
            plt.xlabel('Energy Channel')
            plt.ylabel('Counts/bin')
            plt.yscale('log')
            plt.hist(energies, bins=energy_bins[0:bin_plot_edge], color='r', rwidth=0.5, zorder=1)

            # Plots the energy bins corresponding to the desired energies as vertical lines
            if flagged_indices.size > 0:
                plt.vlines(energy_bins[flagged_indices.astype(int)], 0, 1e6, zorder=2, alpha=0.75)

            # Saves the figure
            sp_path = f'{sm.results_loc()}Results/{self.unit}/{self.full_day_string}/'
            sm.path_maker(sp_path)
            plt.savefig(f'{sp_path}{scintillator}_Spectrum.png', dpi=500)
            plt.clf()

        if self.template:
            exit()

        return lp_energies, nai_energies

    def data_importer(self):
        """Imports data from data files into arrays and then updates them into the nested dictionary structure."""
        for i in self.scintillators:
            scintillator = self.scintillators[i]
            eRC = scintillator['eRC']
            sm.print_logger('\n', self.log)
            sm.print_logger(f'For eRC {eRC} ({i}):', self.log)
            # Here in case the data files in a custom location are grouped into daily folders
            try:
                if self.THOR or self.SANTIS:
                    complete_filelist = glob.glob(f'{self.import_path}/eRC{eRC}*_lm_{self.full_day_string}_*')
                else:
                    complete_filelist = glob.glob(f'{self.import_path}/eRC{eRC}_lm*_{self.full_day_string}_*')

                assert len(complete_filelist) > 0, 'Empty filelist'

            except AssertionError:
                if self.THOR or self.SANTIS:
                    complete_filelist = glob.glob(f'{self.import_path}/{self.full_day_string}'
                                                  f'/eRC{eRC}*_lm_{self.full_day_string}_*')
                else:
                    complete_filelist = glob.glob(f'{self.import_path}/{self.full_day_string}'
                                                  f'/eRC{eRC}_lm*_{self.full_day_string}_*')

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

            try:
                # Tests to make sure that filelist isn't empty
                assert len(filelist) > 0, 'No files for this scintillator today'

                # Starts actually importing the data
                energy_list = []
                time_list = []

                filetimes = np.array([])
                file_time_gaps = np.array([])
                last_second = 0.0
                files_imported = 0

                print('File|File Behavior|File Time Gap (sec)', file=self.log)
                for day_files in filelist:
                    # Try-except block to log files where GPS and wallclock disagree significantly
                    file_behavior = 'Normal'
                    try:
                        if self.processed:
                            e, t = np.loadtxt(day_files, skiprows=1, usecols=(0, 2), unpack=True)
                            energy_list.append(e)
                            time_list.append(t)
                            filetimes = t
                        else:
                            data = dr.fileNameToData(day_files)
                            if 'energies' in data.columns:
                                energy_list.append(data['energies'].to_numpy())
                            else:
                                energy_list.append(data['energy'].to_numpy())

                            time_list.append(data['SecondsOfDay'].to_numpy())
                            filetimes = data['SecondsOfDay'].to_numpy()

                    except Exception as ex:
                        if str(ex) == 'wallclock and GPS clocks in significant disagreement':
                            file_behavior = 'Disagreement'
                            pass
                        else:  # Mostly here so that if the reader ever runs into other errors I'll know about them
                            print(f'line {ex.__traceback__.tb_lineno}:')
                            print(f'{type(ex).__name__}: {ex}')
                            exit()

                    # Determines the time gaps between adjacent files
                    first_second = filetimes[0]
                    file_time_gap = first_second - last_second if files_imported > 0 else 0
                    file_time_gaps = np.append(file_time_gaps, file_time_gap)
                    last_second = filetimes[-1]
                    files_imported += 1

                    print(f'{files_imported}/{len(filelist)} files imported', end='\r')
                    print(f'{day_files}|{file_behavior}|{file_time_gap}', file=self.log)

                # Makes the final arrays and exports them
                scintillator.update({'filelist': filelist})
                times = np.concatenate(time_list)
                # Corrects for the fact that the first 200-300 seconds of the next day are included in the last file
                day_change_array = np.array(np.where(np.diff(times) < -86400))
                if day_change_array.size > 0:
                    change_index = int(day_change_array[0]) + 1
                    times = np.append(times[0:change_index], times[change_index:] + 86400.0)

                scintillator.update({'time': times})
                scintillator.update({'energy': np.concatenate(energy_list)})
                self.scintillators.update({i: scintillator})

                print('\n', file=self.log)
                print(f'Total Counts: {len(np.concatenate(time_list))}', file=self.log)
                print(f'Average time gap: {np.sum(file_time_gaps) / len(file_time_gaps)}', file=self.log)
                print('\n', file=self.log)
            except (AssertionError, ValueError):
                print('Missing data for the specified day.', file=self.log)
                print('Missing data for the specified day.', end='\r')
                continue


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

    def scatterplot_maker(self, timescales, detector, event_number, filelist):
        """Makes the short event scatter plots.

        Parameters
        ----------
        timescales : list
            A list of the timescales (in seconds) that the scatter plots are generated in.
        detector : Detector object
            The detector object used to store all the data and relevant information.
        event_number : int
            A number corresponding to the event's number for that particular scintillator (i.e. 1 would be the first
            event for whatever scintillator, 2 would be the second, and so on).
        filelist : list
            A list of all the files to be searched when looking for the beginning of an event. Once the name of the
            file for an event has been found it is added to the scatter plot title.

        Returns
        -------
        list
            Returns a shorter version of the filelist for the requested scintillator. This new list is
            eventually used when the next scatter plot is generated to make file finding faster.

        """
        times = detector.attribute_retriever(self.scintillator, 'time')
        energies = detector.attribute_retriever(self.scintillator, 'energy')

        event_times = times[self.start:self.stop]
        event_energies = energies[self.start:self.stop]
        event_length = event_times[-1] - event_times[0]
        event_time = times[self.start] - 86400 if times[self.start] > 86400 else times[self.start]

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
            ax.scatter(times, energies + 0.6, s=dot_size, zorder=1, alpha=1.0)
            ax.set_xlabel(f'Time (Seconds, {ts}s total)')
            ax.set_ylabel('Energy Channel')
            # Lines appear (100*percent)% to the left or right of event start/stop depending on subplot timescale
            percent = 0.001
            ax.vlines([event_times[0] - percent*ts, event_times[-1] + percent*ts], 0, 1e5,
                      colors='r', linewidth=1, zorder=-1, alpha=0.3)

        # Adds the name of the relevant data file to the scatter plot
        event_file = ''
        eventfound = False
        files_examined = 0
        new_filelist = []
        for day_files in filelist:
            if detector.processed:
                e, temp_t_array = np.loadtxt(day_files, skiprows=1, usecols=(0, 2), unpack=True)
            else:
                try:
                    data = dr.fileNameToData(day_files)
                except Exception as ex:
                    if str(ex) == 'wallclock and GPS clocks in significant disagreement':
                        pass
                    else:  # Mostly here so that if the reader ever runs into other errors I'll know about them
                        print(f'line {ex.__traceback__.tb_lineno}:')
                        print(f'{type(ex).__name__}: {ex}')
                        exit()

                temp_t_array = data['SecondsOfDay'].to_numpy()

            for j in temp_t_array:
                if j == event_time:
                    event_file = day_files
                    eventfound = True
                    break

            if eventfound:
                new_filelist = filelist[files_examined:]
                break

            files_examined += 1
        plt.title(f'Obtained from {event_file}', fontsize=15, y=-0.4)

        # Saves the scatter plot
        # Note: with this code, if an event happens in that 200-300 seconds of the next day that are included in the
        # last file, the image will have the wrong date in its name (though the timestamp in the title will
        # always be correct)
        scatterpath = f'{sm.results_loc()}Results/{detector.unit}/{detector.full_day_string}/scatterplots/'
        sm.path_maker(scatterpath)
        figure1.savefig(f'{scatterpath}{detector.full_day_string}_{self.scintillator}_event{event_number}.png')
        plt.close(figure1)

        # Makes a json file for the event
        eventpath = f'{sm.results_loc()}Results/{detector.unit}/{detector.full_day_string}/event files/short events/'
        sm.path_maker(eventpath)
        event_frame = pd.DataFrame()
        event_frame['SecondsOfDay'] = event_times
        event_frame['energies'] = event_energies
        event_frame['file'] = event_file  # Note: this column will be filled by the same file name over and over again
        event_frame.to_json(f'{eventpath}{detector.full_day_string}_{self.scintillator}_event{event_number}.json')

        return new_filelist


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

    def highest_zscore(self, z_scores):
        """ Identifies the highest z-score and its corresponding bin for an event."""
        glow_scores = z_scores[self.start:self.stop]
        highest_score = np.max(glow_scores)
        self.peak_index = np.argmax(glow_scores) + self.start
        return highest_score

    def beginning_and_end_seconds(self, bins10sec):
        """ Retrieves the beginning and total length of an event in seconds."""
        glow_times = bins10sec[self.start:self.stop]
        first_sec = glow_times[0]
        length = self.length * 10
        last_sec = first_sec + length
        return first_sec, last_sec
