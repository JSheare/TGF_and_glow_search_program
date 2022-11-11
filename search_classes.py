import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import search_module as sm
import DataReaderFinal as dr

# Classes used in search.py


# Class used to store all relevant information about the detector and the data
class Detector:
    def __init__(self, unit, day, modes):
        self.unit = unit
        self.day = day
        self.modes = modes

        self.THOR = False
        self.GODOT = False
        self.SANTIS = False

        self.custom = False
        self.processed = False
        self.plastics = False
        self.import_path = ''

        if self.unit == 'GODOT':
            self.GODOT = True
            if 'processed' in self.modes:
                self.processed = True
                self.import_path = f'{sm.G_processed_data_loc()}{day[0:4]}/'
            else:
                self.import_path = f'{sm.G_raw_data_loc()}{day}/'

            self.scintillators = {'NaI': {'eRC': '1490', 'filelist': [],
                                          'time': np.array([]), 'energy': np.array([])},
                                  'LP': {'eRC': '1491', 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}

        elif self.unit[0:4] == 'THOR':
            self.THOR = True
            self.import_path = f'{sm.T_raw_data_loc()}{unit}/Data/{day}/'
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
            self.import_path = f'{sm.S_raw_data_loc()}{day}/'
            self.scintillators = {'LP': {'eRC': '2549', 'filelist': [],
                                         'time': np.array([]), 'energy': np.array([])}}
        else:
            print('Not a valid detector')
            exit()

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

    # Retrieves the requested attribute for a particular scintillator
    def attribute_retriever(self, scintillator, attribute):  # Make it possible to request multiple attributes at once?
        medium = self.scintillators[scintillator]
        desired_attribute = medium[attribute]
        return desired_attribute

    # Combines the time and energy arrays for multiple scintillators
    def scintillator_combiner(self, *args):
        time = np.array([])
        energy = np.array([])
        iterations = 0
        for arg in args:
            new_time = self.attribute_retriever(arg, 'time')
            new_energy = self.attribute_retriever(arg, 'energy')
            if iterations == 0:
                time = new_time
                energy = new_energy
            else:
                time = np.append(time, new_time)
                sorting_order = np.argsort(time)
                time = np.sort(time)
                energy = np.append(energy, new_energy)
                energy = energy[sorting_order.astype(int)]

        return time, energy

    # Imports data from datafiles into arrays
    def data_importer(self, datetime_logs):
        for i in self.scintillators:
            scintillator = self.scintillators[i]
            eRC = scintillator['eRC']
            sm.print_logger('\n', datetime_logs)
            sm.print_logger(f'For eRC {eRC} ({i}):', datetime_logs)
            if self.THOR:
                complete_filelist = glob.glob(f'{self.import_path}eRC{eRC}*_lm_{self.day}_*')
            else:
                complete_filelist = glob.glob(f'{self.import_path}eRC{eRC}_lm4_*_{self.day}_*')

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
                scintillator.update({'time': np.concatenate(time_list)})
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
    def scatterplot_maker(self, timescales, filelist, times, energies, event_number, date_timestamp, unit, mode):
        # so many arguments
        event_times = times[self.start:self.stop]
        event_energies = energies[self.start:self.stop]
        event_length = event_times[-1] - event_times[0]

        figure1 = plt.figure(figsize=[20, 11.0])
        figure1.suptitle(f'{self.scintillator} Event {str(event_number)} for {date_timestamp}, {len(event_energies)} '
                         f'counts', fontsize=20)
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
        event_time = times[self.start]
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
        full_day_string = f'{str(date_timestamp)[2:4]}{str(date_timestamp)[5:7]}{str(date_timestamp)[8:10]}'  # Lazy fix
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
