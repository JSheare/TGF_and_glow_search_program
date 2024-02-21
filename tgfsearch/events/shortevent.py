"""A class used to keep track of and process short events."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

import tgfsearch.tools as tl
import tgfsearch.parameters as params


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
    len_subscore : int
        The event's length score given by the ranking system.
    clumpiness_subscore : int
        The event's clumpiness score given by the ranking system.
    hel_subscore : int
        The event's high energy leading count score given by the ranking system.
    weather_subscore : int
        The event's weather score given by the ranking system.
    total_score : int
        The event's total score as determined by the subscores. Calculated by the ranking system.
    rank : int
        The event's rank among all the other events found.

    """

    def __init__(self, event_start, event_length, scintillator):
        self.start = int(event_start)
        self.length = int(event_length)
        self.stop = int(event_start + event_length)
        self.scintillator = scintillator

        # Ranking scores
        self.len_subscore = 0
        self.clumpiness_subscore = 0
        self.hel_subscore = 0
        self.weather_subscore = 0
        self.total_score = 0

        self.rank = 0

    # String casting overload
    def __str__(self):
        return f'{self.scintillator} short event at index {self.start}'

    # Debugging string dunder
    def __repr__(self):
        return f'{self.scintillator} short event; start:{self.start}; stop:{self.stop}; length:{self.length}'

    def _calculate_subscores(self, weather_cache, detector, times, energies):
        """Calculates the event's ranking subscores and updates them."""
        clumpiness = 0
        clump_counts = 0

        high_energy_lead = 0
        leading_counts = 0

        for i in range(self.start, self.stop):
            if i > self.start:
                difference = times[i] - times[i - 1]
                if difference < params.DIFFERENCE_THRESH:
                    clump_counts += 1
                    if clump_counts == 1:
                        leading_counts += 1
                        if energies[i - 1] >= params.HIGH_ENERGY_LEAD_THRESH:  # Clump starts on the previous index
                            high_energy_lead += 1

                    # For the last count in the event
                    if i == self.stop - 1:
                        clump_counts += 1

                else:
                    # Adding to clumpiness when there's a clump of three or more
                    if clump_counts >= 3:
                        clumpiness += 1
                        clump_counts = 0

                    # Adding to the clumpiness when the gap between sufficient counts is greater than the threshold
                    if difference >= params.GAP_THRESH and clump_counts == 0:
                        clumpiness += 1

                    clump_counts = 0

        clumpiness /= self.length
        if leading_counts > 0:
            high_energy_lead /= leading_counts

        # Calculating the length subscore
        self.len_subscore = 1 / params.GOOD_LEN_THRESH * self.length

        # Calculating the clumpiness subscore
        self.clumpiness_subscore = 1 / (1 + np.e ** (40 * (clumpiness - params.CLUMPINESS_TOSSUP)))

        # Calculating high energy leading count subscore
        self.hel_subscore = np.e ** -(5.3 * high_energy_lead)

        # Getting weather subscore
        if detector.location['Nearest weather station'] != '':
            local_date, local_time = tl.convert_to_local(detector.full_date_str, times[self.start])
            self.weather_subscore = tl.get_weather_conditions(local_date, local_time, detector, weather_cache)
        else:
            self.weather_subscore = -1

    def calculate_score(self, weather_cache, detector, times, energies):
        """Calculates the event's score (likelihood of being an interesting event) and updates its total_score
        attribute.

        Parameters
        ----------
        weather_cache : dict
            A cache containing weather info for the day being investigated. If no info has been added yet, this
            function will find and add it during the scoring process.
        detector : Detector object
            The detector object used to store all data and relevant information
        times : np.array
            A numpy array containing times for each count.
        energies : np.array
            A numpy array containingenergies for each count.

        """

        self._calculate_subscores(weather_cache, detector, times, energies)
        # If weather info couldn't be obtained, the weather subscore is removed and the remaining weights are adjusted
        # so that they stay proportional to one another
        if self.weather_subscore == -1:
            proportionality = 1 / (params.LEN_WEIGHT + params.CLUMPINESS_WEIGHT + params.HEL_WEIGHT)
            self.total_score = proportionality * (params.LEN_WEIGHT * self.len_subscore +
                                                  params.CLUMPINESS_WEIGHT * self.clumpiness_subscore +
                                                  params.HEL_WEIGHT * self.hel_subscore)
        else:
            self.total_score = (params.LEN_WEIGHT * self.len_subscore + 
                                params.CLUMPINESS_WEIGHT * self.clumpiness_subscore +
                                params.HEL_WEIGHT * self.hel_subscore +
                                params.WEATHER_WEIGHT * self.weather_subscore)

    def get_filenames(self, filelist_dict, filetime_extrema_dict, times, count_scints=None):
        """Gets the names of the files that the event occurred in.

        Parameters
        ----------
        filelist_dict : dict
            A dictionary containing a list of files for each scintillator that contributed to the event.
        filetime_extrema_dict : dict
            A dictionary containing a list of file time extrema for each scintillator that contributed to the event.
        times : np.array
            A numpy array containing times for each count.
        count_scints : np.array
            A numpy array containing the scintillator each count is from.

        Returns
        -------
        dict
            event_file_dict
                A dictionary containing the name of the files that the event occurred in for each scintillator.
                Normally, this will only contain one entry, but in combo mode there will be an entry for each
                scintillator that contributed to the event.
            filelist_dict
                An updated version of the parameter filelist_dict that only contains files which haven't been
                searched and the most recently searched file.
            filetime_extrema_dict
                An updated version of the parameter filetime_extrema_dict that only contains extrema for files
                which haven't been searched and the most recently searched file.

        """

        event_file_dict = dict()
        for i in range(self.start, self.stop):
            if count_scints is not None:
                scintillator = count_scints[i]
            else:
                scintillator = self.scintillator

            if scintillator not in event_file_dict:
                event_time = times[i] - 86400 if times[i] > 86400 else times[i]
                event_file = ''
                files_examined = 0
                filetime_extrema = filetime_extrema_dict[scintillator]
                filelist = filelist_dict[scintillator]
                for j in range(len(filetime_extrema)):
                    first_time = filetime_extrema[j][0]
                    last_time = filetime_extrema[j][1]
                    if first_time <= event_time <= last_time:
                        event_file = filelist[j]
                        filelist_dict[scintillator] = filelist[files_examined:]
                        filetime_extrema_dict[scintillator] = filetime_extrema[files_examined:]
                        break

                    files_examined += 1

                event_file_dict[scintillator] = event_file

            # So that we don't loop through the whole event for no reason when not in combo mode
            if count_scints is None:
                break

        return event_file_dict, filelist_dict, filetime_extrema_dict

    def make_json(self, event_number, event_file_dict, detector, times, energies, wallclock, count_scints=None):
        """Makes a json file for the event.

        Parameters
        ----------
        event_number : int
            A number corresponding to the event's number for that particular scintillator (i.e. 1 would be the first
            event for whatever scintillator, 2 would be the second, and so on).
        event_file_dict : dict
            A dictionary containing the names of the files that the event occurred in for each scintillator.
            Normally, this will only contain one entry, but in combo mode there will be an entry for each
            scintillator that contributed to the event.
        detector : Detector object
            The detector object used to store all the data and relevant information.
        times : np.array
            A numpy array containing times for each count.
        energies : np.array
            A numpy array containing energies for each count.
        wallclock : np.array
            A numpy array containing wallclock times for each count.
        count_scints : np.array / None
            A numpy array containing the scintillator each count is from. Should be None unless in combo mode.

        """

        eventpath = (f'{detector.get_results_loc()}/Results/{detector.unit}/'
                     f'{detector.date_str}/event files/short events/')
        tl.make_path(eventpath)
        event_frame = pd.DataFrame()
        event_frame['wc'] = wallclock[self.start:self.stop]
        event_frame['SecondsOfDay'] = times[self.start:self.stop]
        event_frame['energies'] = energies[self.start:self.stop]
        event_frame['count_scintillator'] = count_scints[self.start:self.stop] if count_scints is not None else (
                [self.scintillator] * self.length)
        # Note: the column will be filled by the same string over and over again
        event_frame['file'] = ', '.join([event_file_dict[scintillator] for scintillator in event_file_dict])

        # Saves the json file
        event_num_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(event_number)))
        rank_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(self.rank)))
        event_frame.to_json(f'{eventpath}{detector.date_str}_{"CM" if count_scints is not None else self.scintillator}_'
                            f'event{event_num_padding}{event_number}_rank{rank_padding}{self.rank}.json')

    def make_scatterplot(self, event_number, event_file_dict, detector, times, energies, count_scints=None):
        """Makes the short event scatter plots.

        Parameters
        ----------
        event_number : int
            A number corresponding to the event's number for that particular scintillator (i.e. 1 would be the first
            event for whatever scintillator, 2 would be the second, and so on).
        event_file_dict : dict
            A dictionary containing the names of the files that the event occurred in for each scintillator.
            Normally, this will only contain one entry, but in combo mode there will be an entry for each
            scintillator that contributed to the event.
        detector : Detector object
            The detector object used to store all the data and relevant information.
        times : np.array
            A numpy array containing times for each count.
        energies : np.array
            A numpy array containing energies for each count.
        count_scints : np.array / None
            A numpy array containing the scintillator each count is from. Should be None unless in combo mode.

        """

        # Subplot timescales
        timescales = [1e-4, 0.005, 2]  # 100 microseconds, 5 milliseconds, 2 seconds

        # Dot colors. Note that if an instrument with more than just NaI, SP, MP, and LP is ever added, this
        # will result in key errors
        colors = {'NaI': 'b', 'SP': 'm', 'MP': 'g', 'LP': 'darkgoldenrod'}

        # Truncated time and energy arrays to speed up scatter plot making
        fraction_of_day = 1 / 64
        spacer = int((len(times) * fraction_of_day) / 2)
        left_edge = 0 if self.start - spacer < 0 else self.start - spacer
        right_edge = (len(times) - 1) if self.stop + spacer > (len(times) - 1) else self.stop + spacer
        if count_scints is not None:
            times_dict, energies_dict = tl.separate_data(times, energies, count_scints, left_edge, right_edge)
        else:
            times_dict = {self.scintillator: times[left_edge:right_edge]}
            energies_dict = {self.scintillator: energies[left_edge:right_edge]}

        event_times = times[self.start:self.stop]
        event_length = event_times[-1] - event_times[0]  # in seconds

        figure1 = plt.figure(figsize=[20, 11.0], dpi=150.)
        figure1.suptitle(f'{"CM" if count_scints is not None else self.scintillator} Event {str(event_number)}, '
                         f'{dt.datetime.utcfromtimestamp(times[self.start] + detector.first_sec)} UTC, '
                         f'{self.length} counts \n Weather: {tl.weather_from_score(self.weather_subscore)} \n'
                         f'Rank: {self.rank}', fontsize=20)
        ax1 = figure1.add_subplot(3, 1, 1)
        ax2 = figure1.add_subplot(3, 1, 2)
        ax3 = figure1.add_subplot(3, 1, 3)
        ax_list = [ax1, ax2, ax3]
        assert len(ax_list) == len(timescales)

        for i in range(len(ax_list)):
            ts = timescales[i]
            ax = ax_list[i]
            padding = (ts - event_length) / 2
            if event_length >= ts:
                best_time = event_times[np.argmin(np.abs(event_times - np.roll(event_times, 1)))]
                ax.set_xlim(xmin=best_time - (ts / 2), xmax=best_time + (ts / 2))
            else:
                ax.set_xlim(xmin=event_times[0] - padding, xmax=event_times[-1] + padding)

            dot_size = 5 if ts == timescales[0] else 3  # makes larger dots for top plot
            ax.set_yscale('log')
            ax.set_ylim([0.5, 1e5])
            for scintillator in times_dict:
                ax.scatter(times_dict[scintillator], energies_dict[scintillator] + 0.6,
                           s=dot_size, zorder=1, alpha=0.5, label=scintillator, color=colors[scintillator])

            ax.set_xlabel(f'Time (Seconds, {ts}s total)')
            ax.set_ylabel('Energy Channel')
            # Lines appear (100*percent)% to the left or right of event start/stop depending on subplot timescale
            percent = 0.001
            # Event start line is orange, end line is red
            ax.vlines([event_times[0] - percent*ts, event_times[-1] + percent*ts], 0, 1e5,
                      colors=['orange', 'r'], linewidth=1, zorder=-1, alpha=0.3)

        # Adds a legend to the plot if we're in combo mode
        if count_scints is not None:
            plt.legend(loc='lower right')

        # Adds the name of the relevant data file to the scatter plot
        if count_scints is None:
            plt.title(f'Obtained from {event_file_dict[self.scintillator]}', fontsize=15, y=-0.4)

        # Saves the scatter plot
        # Note: with this code, if an event happens in that 200-300 seconds of the next day that are included in the
        # last file, the image will have the wrong date in its name (though the timestamp in the scatter plot title will
        # always be correct)
        scatterpath = (f'{detector.get_results_loc()}/Results/{detector.unit}/'
                       f'{detector.date_str}/scatterplots/')
        tl.make_path(scatterpath)
        event_num_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(event_number)))
        rank_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(self.rank)))
        figure1.savefig(f'{scatterpath}{detector.date_str}_{"CM" if count_scints is not None else self.scintillator}_'
                        f'event{event_num_padding}{event_number}_rank{rank_padding}{self.rank}.png')
        plt.close(figure1)
