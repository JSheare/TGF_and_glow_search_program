import glob
import numpy as np
import pandas as pd
import scipy.signal as signal
import matplotlib.pyplot as plt
import search_module as sm
import DataReaderFinal as dr
import datetime
from datetime import timedelta as td

# Classes used in search.py


# Class used to store all relevant information about the detector and the data
class Detector:
    def __init__(self, unit, day, modes):
        self.unit = unit
        self.day = day
        self.location = 'location'  # Eventually this will be fetched from a function in search_module
        self.modes = modes

        self.THOR = False
        self.GODOT = False
        self.SANTIS = False

        self.custom = False
        self.processed = False
        self.plastics = False
        self.template = False

        self.import_path = ''
        self.good_lp_calibration = False

        if self.unit == 'GODOT':
            self.GODOT = True
            if 'processed' in self.modes:
                self.processed = True
                self.import_path = f'{sm.G_processed_data_loc()}/{day[0:4]}/'
            else:
                self.import_path = f'{sm.G_raw_data_loc()}/{day}/'

            self.scintillators = {'NaI': {'eRC': '1490', 'filelist': [],
                                          'time': np.array([]), 'energy': np.array([])},
                                  'LP': {'eRC': '1491', 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}

        elif self.unit[0:4] == 'THOR':
            self.THOR = True
            self.import_path = f'{sm.T_raw_data_loc()}/{unit}/Data/{day}/'
            self.scintillators = {'NaI': {'eRC': sm.T_eRC(self.unit)[0], 'filelist': [],
                                          'time': np.array([]), 'energy': np.array([])},
                                  'SP': {'eRC': sm.T_eRC(self.unit)[1], 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])},
                                  'MP': {'eRC': sm.T_eRC(self.unit)[2], 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])},
                                  'LP': {'eRC': sm.T_eRC(self.unit)[3], 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}

        elif self.unit == 'SANTIS':
            self.SANTIS = True
            self.import_path = f'{sm.S_raw_data_loc()}/{day}/'
            self.scintillators = {'LP': {'eRC': '2549', 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}
        else:
            print('Not a valid detector')
            exit()

        if 'custom' in self.modes:
            self.custom = True  # Not necessary for anything right now
            self.import_path = f'{sm.C_raw_data_loc()}/'

        if 'processed' in self.modes:
            self.processed = True
            if self.unit != 'GODOT':
                print('Processed data is only accessible for GODOT')
                exit()

        if 'plastics' in self.modes:
            self.plastics = True

        if 'template' in self.modes:
            self.template = True

    # Retrieves the requested attribute for a particular scintillator
    def attribute_retriever(self, scintillator, attribute):  # Make it possible to request multiple attributes at once?
        medium = self.scintillators[scintillator]
        desired_attribute = medium[attribute]
        return desired_attribute

    # Makes the energy spectra histograms for the LP and NaI scintillators
    def spectra_maker(self, date_timestamp, full_day_string, log):
        lp_energies = np.array([])
        nai_energies = np.array([])
        if self.THOR:
            bin_range = 65535.0
            bin_number = 65536
            band_starts = [1800, 5600]
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
                sm.print_logger('Entering template mode...', log)
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
                except FileNotFoundError:  # Here just in case there is no template found for the specified location
                    # This LP calibration is likely to be pretty inaccurate
                    sm.print_logger('No LP template found for this location', log)

                lp_energies = energy_bins[flagged_indices.astype(int)]

            else:  # NaI scintillator
                # Takes the sum of each bin with its two closest neighboring bins on either side
                sums = energy_hist
                for i in range(2):
                    sums += (np.roll(energy_hist, i + 1) + np.roll(energy_hist, i - 1))

                # Looks for the location of the maximum sum within the two bands where the peaks are likely to be
                for th in range(len(band_starts)):
                    band_max = np.argmax(sums[band_starts[th]:band_ends[th]]) + int(band_starts[th])
                    flagged_indices = np.append(flagged_indices, band_max)
                nai_energies = energy_bins[flagged_indices.astype(int)]

            # Plots the actual spectrum
            plt.figure(figsize=[20, 11.0])
            plt.title(f'Energy Spectrum for {scintillator}, {str(date_timestamp)}', loc='center')
            plt.xlabel('Energy Channel')
            plt.ylabel('Counts/bin')
            plt.yscale('log')
            plt.hist(energies, bins=energy_bins[0:bin_plot_edge], color='r', rwidth=0.5, zorder=1)

            # Saves energy bins corresponding to the desired energies and plots them as vertical lines
            plt.vlines(energy_bins[flagged_indices.astype(int)], 0, 1e6, zorder=2, alpha=0.75)

            # Saves the figure
            sp_path = f'{sm.results_loc()}Results/{self.unit}/{full_day_string}/'
            sm.path_maker(sp_path)
            plt.savefig(f'{sp_path}{scintillator}_Spectrum.png', dpi=500)
            plt.clf()

        if self.template:
            exit()

        return lp_energies, nai_energies

    # Imports data from datafiles into arrays
    def data_importer(self, datetime_logs):
        for i in self.scintillators:
            scintillator = self.scintillators[i]
            eRC = scintillator['eRC']
            sm.print_logger('\n', datetime_logs)
            sm.print_logger(f'For eRC {eRC} ({i}):', datetime_logs)
            # Here in case the data files in a custom location are grouped into daily folders
            try:
                if self.THOR:
                    complete_filelist = glob.glob(f'{self.import_path}/eRC{eRC}*_lm_{self.day}_*')
                else:
                    complete_filelist = glob.glob(f'{self.import_path}/eRC{eRC}_lm4_*_{self.day}_*')

                assert len(complete_filelist) > 0, 'Empty filelist'

            except AssertionError:
                if self.THOR:
                    complete_filelist = glob.glob(f'{self.import_path}{self.day}/eRC{eRC}*_lm_{self.day}_*')
                else:
                    complete_filelist = glob.glob(f'{self.import_path}{self.day}/eRC{eRC}_lm4_*_{self.day}_*')

            # Filters out trace mode files and .txtp files (whatever those are)
            filtered_filelist = []
            unfiltered_filelist = []
            for u in range(len(complete_filelist)):
                current_file = complete_filelist[u]
                if current_file[-3:] == 'xtr' or current_file[-4:] == 'txtp' or current_file[-5:] == 'xtrpp':
                    continue
                elif current_file[-7:] == '.txt.gz':
                    filtered_filelist.append(current_file)
                else:
                    unfiltered_filelist.append(current_file)

            # Eliminates duplicate files
            for h in range(len(unfiltered_filelist)):
                current_file = unfiltered_filelist[h]
                if f'{current_file}.gz' in filtered_filelist:
                    continue
                else:
                    filtered_filelist.append(current_file)

            # Puts the files in the correct order
            files = []
            extensions = []
            for ik in range(len(filtered_filelist)):
                file = filtered_filelist[ik]
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
            for it in range(len(files)):
                filelist.append(f'{files[it]}{extensions[it]}')

            try:
                # Tests to make sure that new_filelist isn't empty
                assert len(filelist) > 0, 'No files for this scintillator today'

                # Starts actually importing the data
                energy_list = []
                time_list = []

                filetimes = np.array([])
                file_time_gaps = np.array([])
                last_second = 0.0
                files_imported = 0

                print('File|File Behavior|File Time Gap (sec)', file=datetime_logs)
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
                    print(f'{day_files}|{file_behavior}|{file_time_gap}', file=datetime_logs)

                # Makes the final arrays and exports them
                scintillator.update({'filelist': filelist})
                times = np.concatenate(time_list)
                # Corrects for the fact that the first 200-300 seconds of the next day are included in the last file
                day_change_array = np.array(np.where(np.diff(times) < -80000))
                if day_change_array.size > 0:
                    change_index = int(day_change_array[0]) + 1
                    times = np.append(times[0:change_index], times[change_index:] + 86400.0)

                scintillator.update({'time': times})
                scintillator.update({'energy': np.concatenate(energy_list)})
                self.scintillators.update({i: scintillator})

                print('\n', file=datetime_logs)
                print(f'Total Counts: {len(np.concatenate(time_list))}', file=datetime_logs)
                print(f'Average time gap: {np.sum(file_time_gaps) / len(file_time_gaps)}', file=datetime_logs)
                print('\n', file=datetime_logs)
            except (AssertionError, ValueError):
                print('Missing data for the specified day.', file=datetime_logs)
                print('Missing data for the specified day.', end='\r')
                continue


class ShortEvent:
    def __init__(self, event_start, event_length, scintillator):
        self.start = int(event_start)
        self.length = int(event_length)
        self.stop = int(event_start + event_length)
        self.scintillator = scintillator

    # Makes the short event scatter plots
    def scatterplot_maker(self, timescales, filelist, times, energies, event_number, first_sec, timestamp, unit, mode):
        # so many arguments
        # I should probably just make some of this stuff attributes of the detector object
        event_times = times[self.start:self.stop]
        event_energies = energies[self.start:self.stop]
        event_length = event_times[-1] - event_times[0]
        event_time = times[self.start] - 86400 if times[self.start] > 86400 else times[self.start]

        figure1 = plt.figure(figsize=[20, 11.0])
        figure1.suptitle(f'{self.scintillator} Event {str(event_number)}, '
                         f'{datetime.datetime.utcfromtimestamp(times[self.start] + first_sec)} UTC, '
                         f'{len(event_energies)} counts', fontsize=20)
        ax1 = figure1.add_subplot(3, 1, 1)
        ax2 = figure1.add_subplot(3, 1, 2)
        ax3 = figure1.add_subplot(3, 1, 3)
        ax_list = [ax1, ax2, ax3]

        for v in range(len(ax_list)):
            ts = timescales[v]
            ax = ax_list[v]
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
            if mode == 'processed':
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
        # last file, the image file will have the wrong date in its name. Though the timestamp in the title will
        # always be correct.
        full_day_string = f'{str(timestamp)[2:4]}{str(timestamp)[5:7]}{str(timestamp)[8:10]}'  # Lazy fix
        scatterpath = f'{sm.results_loc()}Results/{unit}/{full_day_string}/scatterplots/'
        sm.path_maker(scatterpath)
        figure1.savefig(f'{scatterpath}{full_day_string}_{self.scintillator}_event{event_number}.png')
        plt.close(figure1)
        return new_filelist


class PotentialGlow:
    def __init__(self, glow_start, glow_length):
        self.start = int(glow_start)
        self.length = int(glow_length)
        self.stop = int(glow_start + glow_length - 1) if self.length > 1 else int(glow_start + glow_length)
        self.peak_index = 0

    def highest_zscore(self, z_scores):
        glow_scores = z_scores[self.start:self.stop]
        highest_score = np.max(glow_scores)
        self.peak_index = np.argmax(glow_scores) + self.start
        return highest_score

    def glow_length_and_beginning_seconds(self, bins10sec):
        glow_times = bins10sec[self.start:self.stop]
        first_sec = glow_times[0]
        length = self.length * 10
        return first_sec, length
