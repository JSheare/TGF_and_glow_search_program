import sys
import warnings
import datetime as dt
import pickle as pickle
import glob as glob
import os as os
import psutil as psutil
import gc as gc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
import search_classes as sc
import search_module as sm


def long_event_cutoff(detector_obj, chunk_obj=None):
    operating_obj = chunk_obj if chunk_obj is not None else detector_obj
    # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
    le_times = np.array([])
    le_energies = np.array([])
    for scintillator in operating_obj.long_event_scint_list:
        existing_calibration = True if len(
            detector_obj.attribute_retriever(scintillator, 'calibration')) == 2 else False
        if scintillator == operating_obj.long_event_scint_list[0]:
            le_times = operating_obj.attribute_retriever(scintillator, 'time')
            if skcali or not existing_calibration:
                le_energies = operating_obj.attribute_retriever(scintillator, 'energy')
            else:
                calibration = detector_obj.attribute_retriever(scintillator, 'calibration')
                le_energies = sm.channel_to_mev(operating_obj.attribute_retriever(scintillator, 'energy'), calibration,
                                                scintillator)
        else:
            le_times = np.append(le_times, operating_obj.attribute_retriever(scintillator, 'time'))
            if skcali or not existing_calibration:
                le_energies = np.append(le_energies, operating_obj.attribute_retriever(scintillator, 'energy'))
            else:
                calibration = detector_obj.attribute_retriever(scintillator, 'calibration')
                le_energies = np.append(le_energies, sm.channel_to_mev(operating_obj.attribute_retriever(
                    scintillator, 'energy'), calibration, scintillator))

    # Removes entries that are below a certain cutoff energy
    if (not detector_obj.good_lp_calibration and 'LP' in detector_obj.long_event_scint_list) or skcali:
        print('Missing calibration(s), beware radon washout!',
              file=detector_obj.log)
    else:
        energy_cutoff = 1.9  # MeV
        cut_indices = np.where(le_energies < energy_cutoff)
        le_times = np.delete(le_times, cut_indices)

    if detector.processed:
        le_times = le_times - first_sec

    return le_times


def short_event_search(detector_obj, prev_event_numbers=None, low_mem=False):
    # Parameters:
    rollgap = 18 if aircraft else 4
    event_time_spacing = 1e-3  # 1 millisecond
    event_min_counts = 10
    noise_cutoff_energy = 300
    min_noise_counts = 3
    channel_range_width = 300
    channel_ratio = 0.5

    event_numbers = prev_event_numbers if prev_event_numbers is not None else {}
    for scintillator in detector_obj:
        if scintillator != 'LP' and not allscints:
            continue

        sm.print_logger(f'Searching eRC {detector_obj.attribute_retriever(scintillator, "eRC")} ({scintillator})...',
                        detector_obj.log)
        times = detector_obj.attribute_retriever(scintillator, 'time')
        energies = detector_obj.attribute_retriever(scintillator, 'energy')
        filelist = detector_obj.attribute_retriever(scintillator, 'filelist')
        filetime_extrema = detector_obj.attribute_retriever(scintillator, 'filetime_extrema')

        # Checks for an event by looking for a certain number of counts (rollgap + 1) in a small timeframe
        potential_event_list = []
        event_start = 0
        event_length = 0

        # Stats
        total_potential_events = 0
        total_threshold_reached = 0
        removed_len = 0
        removed_channel_ratio = 0
        removed_low_energy = 0
        removed_noise_general = 0
        removed_crs = 0

        rolled_array = np.roll(times, -rollgap)
        interval = rolled_array - times
        for i in range(len(interval)):
            if 0 < interval[i] <= event_time_spacing:
                # Records the beginning index of a potential event
                if event_length == 0:
                    event_start = i
                    event_length = 1 + rollgap  # 1 for first count, rollgap for the others
                    total_potential_events += 1
                # Measures the length of a potential event
                else:
                    event_length += 1

                # Counts the total number of times that the detection threshold was reached
                total_threshold_reached += 1

            # Records the rough length of a potential event
            if interval[i] > event_time_spacing and event_length > 0:
                # Keeps potential event if it's longer than the specified minimum number of counts
                if event_length >= event_min_counts:
                    potential_event = sc.ShortEvent(event_start, event_length, scintillator)
                    potential_event_list.append(potential_event)

                if event_length < event_min_counts:
                    removed_len += 1

                event_start = 0
                event_length = 0

        # Eliminates noisy events
        f_potential_event_list = []
        for event in potential_event_list:
            event_energies = energies[event.start:event.stop]
            event_length = len(event_energies)

            # Checks that there are fewer counts in the higher energy channels than the low energy channels
            # High/low channel ratio is not checked for THOR
            good_channel_ratio = True if detector_obj.THOR else False
            if event_length >= 30:
                low_channel_counts = 0
                high_channel_counts = 0
                for i in event_energies:
                    if 200 <= i <= (200 + channel_range_width):
                        low_channel_counts += 1

                    if (300 + channel_range_width) <= i <= (300 + 2 * channel_range_width):
                        high_channel_counts += 1

                try:
                    if high_channel_counts / low_channel_counts <= channel_ratio:
                        good_channel_ratio = True

                except ZeroDivisionError:
                    pass
            # High/low channel ratio is also not checked for events with < 30 counts
            else:
                good_channel_ratio = True

            # Eliminates the low energy events, which are likely just noise
            # At least min_noise_counts must be above noise_cutoff_energy
            max_indices = (-event_energies).argsort()[:min_noise_counts]
            max_energies = np.array([])
            for i in max_indices:
                max_energies = np.append(max_energies, event_energies[i])

            is_greater_than_thresh = np.all((max_energies > noise_cutoff_energy))

            # Adds everything that passed the two tests to a separate array
            if is_greater_than_thresh and good_channel_ratio:
                f_potential_event_list.append(event)
            else:
                if not good_channel_ratio and is_greater_than_thresh:
                    removed_channel_ratio += 1
                elif not is_greater_than_thresh and good_channel_ratio:
                    removed_low_energy += 1
                else:
                    removed_noise_general += 1

        potential_event_list = f_potential_event_list

        # Eliminates events that are just successive cosmic ray showers
        if aircraft:
            a_potential_event_list = []
            # Parameters
            difference_threshold = 2e-6
            gap_threshold = 10e-6
            clumpiness_threshold = 0.27
            for event in potential_event_list:
                time_slice = times[event.start:event.stop] - times[event.start]
                if event.length >= 30:
                    a_potential_event_list.append(event)
                    continue
                else:
                    clumpiness = 0
                    clump_counts = 0
                    for i in range(len(time_slice) - 1):
                        difference = time_slice[i+1] - time_slice[i]
                        if difference < difference_threshold:
                            clump_counts += 1
                            if i+1 == len(time_slice):
                                clump_counts += 1
                        else:
                            # Adding to clumpiness when there's a clump of three or more
                            if clump_counts >= 3:
                                clumpiness += 1
                                clump_counts = 0

                            # Adding to clumpiness when the gap between sufficient counts is greater than the threshold
                            if difference >= gap_threshold and clump_counts == 0:
                                clumpiness += 1

                            clump_counts = 0

                    clumpiness /= len(time_slice)
                    if clumpiness < clumpiness_threshold:
                        a_potential_event_list.append(event)
                    else:
                        removed_crs += 1

            potential_event_list = a_potential_event_list

        sm.print_logger(f'\n{len(potential_event_list)} potential events recorded', detector_obj.log)
        if removed_len > 0:
            sm.print_logger(f'{removed_len} events removed due to insufficient length', detector_obj.log)

        if removed_channel_ratio > 0:
            sm.print_logger(f'{removed_channel_ratio} events removed due to noise (bad high-to-low channel ratio)',
                            detector_obj.log)

        if removed_low_energy > 0:
            sm.print_logger(f'{removed_low_energy} events removed due to noise (minimum energy threshold not reached)'
                            , detector_obj.log)

        if removed_noise_general > 0:
            sm.print_logger(f'{removed_noise_general} events removed due to noise (general)', detector_obj.log)

        if removed_crs > 0:
            sm.print_logger(f'{removed_crs} events removed due to successive CRS', detector_obj.log)

        sm.print_logger(f'Detection threshold reached {total_threshold_reached} times', detector_obj.log)
        print('\n', file=detector_obj.log)

        if len(potential_event_list) > 0:
            print('\n')
            sm.print_logger('Generating scatter plots and event files...', detector_obj.log)
            print('\n', file=detector_obj.log)
            print('Potential short events:', file=detector_obj.log)
            for event in potential_event_list:
                start_second = times[event.start] - 86400 if times[event.start] > 86400 else times[event.start]
                print(f'{dt.datetime.utcfromtimestamp(times[event.start] + first_sec)} UTC '
                      f'({start_second} seconds of day)', file=detector_obj.log)

            print('\n', file=detector_obj.log)

            # Makes scatter plots of the resulting potential events
            # Subplot timescales
            ts1 = 1e-4   # 100 microseconds
            ts2 = 0.005  # 5 milliseconds
            ts3 = 2      # 100 milliseconds
            ts_list = [ts1, ts2, ts3]

            if prev_event_numbers is not None:
                plots_already_made = event_numbers[scintillator]
            else:
                plots_already_made = 0

            plots_made = 0
            for i in range(len(potential_event_list)):
                print(f'{plots_made}/{len(potential_event_list)}', end='\r')
                event = potential_event_list[i]
                filelist, filetime_extrema = event.scatterplot_maker(ts_list, detector_obj, i + 1 + plots_already_made,
                                                                     filelist, filetime_extrema)
                plots_made += 1
                print(f'{plots_made}/{len(potential_event_list)}', end='\r')

            event_numbers.update({scintillator: plots_made})
        else:
            print('\n')
            event_numbers.update({scintillator: 0})

    if low_mem:
        return event_numbers


def long_event_search(detector_obj, le_times, existing_hist=None, low_mem=False):
    # Makes one bin for every binsize seconds of the day (plus around 300 seconds more for the next day)
    binsize = 4
    day_bins = np.arange(0, 86700 + binsize, binsize)

    # Creates numerical values for histograms using numpy
    if existing_hist is not None:
        hist_allday = existing_hist
    else:
        hist_allday, bins_allday = np.histogram(le_times, bins=day_bins)
        if low_mem:
            return hist_allday

    flag_threshold = 5

    # Calculates mean and z-scores
    hist_allday_nz = hist_allday[hist_allday.nonzero()]
    mue_val = sum(hist_allday_nz) / len(hist_allday_nz)  # Mean number of counts/bin

    # Removes outliers (3 sigma) and recalculates mue
    abs_zscores = np.abs((hist_allday_nz - mue_val) / np.sqrt(mue_val))
    hist_allday_nz = hist_allday_nz[np.argsort(abs_zscores)]
    hist_allday_nz = hist_allday_nz[::-1]
    for i in range(len(abs_zscores)):
        if abs_zscores[i] < 3:
            hist_allday_nz = hist_allday_nz[i:]
            break

    mue_val = sum(hist_allday_nz) / len(hist_allday_nz)

    mue = np.ones(len(day_bins)-1) * mue_val
    sigma = np.ones(len(day_bins)-1) * np.sqrt(mue_val)

    if aircraft:
        fitting_window = 20
        window = 1
        gap = 5
        num_bins = len(day_bins) - 1
        so_poly = lambda t, a, b, c: a * t ** 2 + b * t + c
        line = lambda t, a, b: a*t + b
        window_index = 0
        while window_index < num_bins - 1:
            fit_curve = so_poly
            l_bins = []
            l_counts = []
            r_bins = []
            r_counts = []
            for i in range(fitting_window):
                # Checking left side bins
                left_bin_index = window_index - (i + 1 + gap)
                if not left_bin_index < 0:
                    if hist_allday[left_bin_index] > 0:
                        left_nz = True
                        right_nz = False if hist_allday[left_bin_index + 1] == 0 else True
                        if not left_bin_index - 1 < 0:
                            if hist_allday[left_bin_index - 1] == 0:
                                left_nz = False

                        if left_nz and right_nz:
                            l_bins.append(day_bins[left_bin_index])
                            l_counts.append(hist_allday[left_bin_index])
                        else:
                            fit_curve = line
                    else:
                        fit_curve = line

                # Checking right side bins
                right_bin_index = window_index + gap + i + window
                if not right_bin_index > (num_bins - 1):
                    if hist_allday[right_bin_index] > 0:
                        left_nz = False if hist_allday[right_bin_index - 1] == 0 else True
                        right_nz = True
                        if not right_bin_index + 1 > (num_bins - 1):
                            if hist_allday[right_bin_index + 1] == 0:
                                right_nz = False

                        if left_nz and right_nz:
                            r_bins.append(day_bins[right_bin_index])
                            r_counts.append(hist_allday[right_bin_index])
                        else:
                            fit_curve = line

                    else:
                        fit_curve = line

            if len(l_bins[::-1] + r_bins) < 3:
                window_index += window
                continue
            else:
                params, pcov = curve_fit(fit_curve, l_bins[::-1] + r_bins, l_counts[::-1] + r_counts)
                window_bins = day_bins[window_index:window_index + window]
                if fit_curve == so_poly:
                    mean_curve = so_poly(window_bins, params[0], params[1], params[2])
                else:
                    mean_curve = line(window_bins, params[0], params[1])

                for i in range(len(mean_curve)):
                    if window_index + i < (num_bins - 1):
                        mue[window_index + i] = mean_curve[i]
                        sigma[window_index + i] = np.sqrt(mean_curve[i])

                window_index += window

    z_scores = np.array([])  # The z-scores themselves
    z_flags = np.array([]).astype(int)
    # Flags only those z-scores > flag_threshold
    for i in range(len(hist_allday)):
        z_score = (hist_allday[i] - mue[i]) / sigma[i]
        z_scores = np.append(z_scores, z_score)
        if z_score >= flag_threshold:
            z_flags = np.append(z_flags, i)

    # Sorts z-flags into actual potential glows
    glow_start = 0
    glow_length = 0
    previous_time = 0
    potential_glow_list = []
    for flag in z_flags:
        if glow_length == 0:  # First zscore only
            glow_start = flag
            glow_length += 1
        elif glow_length > 0 and day_bins[flag] - binsize == previous_time:
            glow_length += 1
        elif (glow_length > 0 and day_bins[flag] - binsize > previous_time) or flag == z_flags[-1]:
            # Makes glow object and fills it out
            glow = sc.PotentialGlow(glow_start, glow_length)
            glow.highest_score = glow.highest_zscore(z_scores)
            glow.start_sec, glow.stop_sec = glow.beginning_and_end_seconds(day_bins, binsize)
            potential_glow_list.append(glow)
            glow_start = flag
            glow_length = 1

        previous_time = day_bins[flag]

    # Rejects events whose bins don't have enough counts (aircraft only)
    if aircraft:
        a_potential_glow_list = []
        minimum_counts = 1000
        for glow in potential_glow_list:
            for bin_counts in hist_allday[glow.start:glow.stop]:
                if bin_counts > minimum_counts:
                    a_potential_glow_list.append(glow)
                    break

        potential_glow_list = a_potential_glow_list

    sm.print_logger('Done.', detector_obj.log)
    if len(potential_glow_list) == 0:
        sm.print_logger('\n', detector_obj.log)
        sm.print_logger(f'There were no potential glows on {date_timestamp}', detector_obj.log)
    else:
        # Logs potential glows and sorts them in descending order depending on their highest z-score
        sm.print_logger('\n', detector_obj.log)
        sm.print_logger('Generating event files...', detector_obj.log)
        highest_scores = []
        print('\n', file=detector_obj.log)
        print('Potential glows:', file=detector_obj.log)

        eventpath = f'{detector_obj.results_loc}Results/{detector_obj.unit}/{detector_obj.full_day_string}/' \
                    f'event files/long events/'
        sm.path_maker(eventpath)
        event_number = 1
        files_made = 0
        for i in range(len(potential_glow_list)):
            print(f'{files_made}/{len(potential_glow_list)}', end='\r')
            glow = potential_glow_list[i]
            highest_score = glow.highest_score
            highest_scores.append(highest_score)
            print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + first_sec)} UTC ({glow.start_sec} '
                  f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                  f'{highest_score}', file=detector_obj.log)

            event_file = open(f'{eventpath}{detector_obj.full_day_string}_event{event_number}_zscore'
                              f'{int(highest_score)}.txt', 'w')
            print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + first_sec)} UTC ({glow.start_sec} '
                  f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                  f'{highest_score}', file=event_file)
            for scintillator in detector_obj.long_event_scint_list:
                print(f'{scintillator}:', file=event_file)
                filelist = detector_obj.attribute_retriever(scintillator, 'filelist')
                filetime_extrema = detector_obj.attribute_retriever(scintillator, 'filetime_extrema')
                files_added = 0
                for j in range(len(filetime_extrema)):
                    first_time = filetime_extrema[j][0]
                    last_time = filetime_extrema[j][1]
                    if first_time <= glow.start_sec <= last_time or \
                            first_time <= glow.stop_sec <= last_time:
                        print(filelist[j], file=event_file)
                        files_added += 1
                    else:
                        if files_added > 0:
                            break

            event_file.close()
            event_number += 1
            files_made += 1
            print(f'{files_made}/{len(potential_glow_list)}', end='\r')

        print('\n', file=detector_obj.log)
        sm.print_logger('Done.', detector_obj.log)

        glow_sorting_order = np.argsort(highest_scores)
        glow_sorting_order = glow_sorting_order[::-1]
        potential_glow_list = [potential_glow_list[s] for s in glow_sorting_order]

        # Establishes detector location based on year and detector name
        location = sm.location(unit, year)

        # Plotting the histograms
        sm.print_logger('\n', detector_obj.log)
        sm.print_logger('Generating Histogram...', detector_obj.log)
        figu = plt.figure(figsize=[20, 11.0])
        plt.title(f'{unit} {location}, {str(date_timestamp)}', loc='center')
        plt.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
        ax1 = figu.add_subplot(5, 1, 1)
        ax2 = figu.add_subplot(5, 1, 2)
        ax3 = figu.add_subplot(5, 1, 3)
        ax4 = figu.add_subplot(5, 1, 4)
        ax5 = figu.add_subplot(5, 1, 5)

        ax_list = [ax2, ax3, ax4, ax5]

        # Plots the histogram for the entire day (just the 1st subplot)
        ax1.bar(day_bins[:-1], hist_allday, alpha=0.5, color='r', width=binsize)
        ax1.set_xlabel('Seconds of Day (UT)')
        ax1.set_ylabel('Counts/bin')
        ax1.plot(day_bins[:-1], mue + flag_threshold*sigma, color='blue', linestyle='dashed')

        # Creates legend
        allday_data = mpatches.Patch(color='r', label='All Energies')
        allday_thresh_sigma = mpatches.Patch(color='blue', label=f'{flag_threshold} Sigma Above All Energies',
                                             linestyle='dashed')
        ax1.legend(handles=[allday_data, allday_thresh_sigma, ], bbox_to_anchor=(1.05, 1), loc=1,
                   borderaxespad=0.)
        ax1.grid(True)

        # Makes the histogram subplots
        for i in range(4):
            try:
                glow = potential_glow_list[i]
                glow.hist_subplotter(ax_list[i], day_bins, hist_allday, mue, sigma, flag_threshold)
            except IndexError:
                continue

        plt.tight_layout()
        plt.plot()

        # Saves the histograms:
        hist_path = f'{detector_obj.results_loc}Results/{unit}/{full_day_string}/'
        sm.path_maker(hist_path)
        plt.savefig(f'{hist_path}{full_day_string}_histogram.png', dpi=500)
        plt.close(figu)


# Disables all warnings in the console because they're mostly just annoying
warnings.filterwarnings('ignore')

# Establish date range from argv 1 and 2
first_date = str(sys.argv[1])
second_date = str(sys.argv[2])

# Establish detector unit from argv 3
unit = str(sys.argv[3])

try:
    modes = sys.argv[4:]
except IndexError:
    modes = []
    pass

# Modes for skipping over certain algorithms (mostly to speed up testing)
skcali = True if 'skcali' in modes else False  # Skip detector calibration
skshort = True if 'skshort' in modes else False  # Skip short event search
skglow = True if 'skglow' in modes else False  # SKip long event search

# Aircraft mode
aircraft = True if 'aircraft' in modes else False

# Pickle mode
picklem = True if 'pickle' in modes else False

# All scintillators mode (all the scintillators will be checked by the short event search algorithm)
allscints = True if 'allscints' in modes else False

# Makes a list of all the dates on the requested range
requested_dates = [first_date]
if first_date != second_date:
    date_int = int(first_date)
    century = '20'
    while True:
        date_int += 1
        date_str = str(date_int)
        # Month rollover
        if int(date_str[4:]) > sm.days_per_month(int(date_str[2:4]), int(century + date_str[0:2])):
            date_int = date_int + 100 - (int(date_str[4:]) - 1)
            date_str = str(date_int)

        # Year rollover
        if int(date_str[2:4]) > 12:
            date_int = (date_int//10000 + 1) * 10000 + 101
            date_str = str(date_int)

        # Note: these rollovers won't work properly for dates outside the 21st century :)

        if date_str == second_date:
            requested_dates.append(date_str)
            break
        else:
            requested_dates.append(date_str)

for date in requested_dates:
    full_day_string = str(date)  # In format yymmdd
    day = int(full_day_string[4:])
    month = int(full_day_string[2:4])
    year = int('20' + full_day_string[0:2])
    date_timestamp = dt.date(year, month, day)  # In format yyyy-mm-dd
    print(f'\n{date_timestamp}:')

    # EPOCH time conversions
    first_sec = (dt.datetime(year, month, day, 0, 0) - dt.datetime(1970, 1, 1)).total_seconds()
    print(f'The first second of {year}-{month}-{day} is: {int(first_sec)}')
    print(f'The last second of {year}-{month}-{day} is: {int(first_sec + 86400)}')

    # Initializes the detector object
    print('Importing data...')
    detector = sc.Detector(unit, first_sec, modes)

    # Logs relevant data files and events in a .txt File
    log_path = f'{detector.results_loc}Results/{unit}/{full_day_string}/'
    sm.path_maker(log_path)
    log = open(f'{log_path}log.txt', 'w')
    log.write(f'The first second of {year}-{month}-{day} is: {int(first_sec)}\n')
    log.write(f'The last second of {year}-{month}-{day} is: {int(first_sec + 86400)}\n')

    # Normal operating mode
    try:
        # Imports the data
        if picklem:
            pickle_path = glob.glob(f'{detector.results_loc}Results/{unit}/{full_day_string}/detector.pickle')
            if len(pickle_path) > 0:
                detector_pickle = open(pickle_path[0], 'rb')
                detector = pickle.load(detector_pickle)
                detector_pickle.close()
                detector.log = log
                # The rest of these are for the modes (which might not necessarily be the same
                # for the serialized detector)
                detector.modes = modes
                detector.template = True if 'template' in modes else False
            else:
                detector.log = log
                detector.data_importer()
                detector.log = None  # serializing open file objects results in errors
                detector.regex = None  # serializing anonymous functions results in errors too
                detector.template = False  # Setting this to false ensures no odd behaviors
                detector_pickle = open(f'{detector.results_loc}Results/{unit}/{full_day_string}/detector.pickle',
                                       'wb')
                pickle.dump(detector, detector_pickle)
                detector_pickle.close()
                detector.log = log
        else:
            detector.log = log
            detector.data_importer()

        # Checks to see if there is actually data for the day
        data_present = True
        if 'LP' in detector.long_event_scint_list:
            necessary_scintillators = detector.long_event_scint_list
        else:
            necessary_scintillators = detector.long_event_scint_list + ['LP']

        for scint in necessary_scintillators:
            if len(detector.attribute_retriever(scint, 'time')) == 0:
                data_present = False
                print('\n\n')
                print('\n', file=detector.log)
                sm.print_logger('No/Missing data for specified day.\n', detector.log)
                print('\n')
                break

        if not data_present:
            raise FileNotFoundError

        print('\n\n')
        print('Done.')

        # Calibrates each scintillator
        if not skcali:
            print('\n')
            sm.print_logger('Calibrating scintillators and generating energy spectra...', detector.log)

            # Calling the calibration algorithm
            detector.spectra_maker()

            sm.print_logger('Done.', detector.log)

        # Short event search
        if not skshort:
            sm.print_logger('\n', detector.log)
            sm.print_logger('Starting search for short events...', detector.log)
            sm.print_logger('\n', detector.log)

            # Calling the short event search algorithm
            short_event_search(detector)

            sm.print_logger('Done.', detector.log)
            sm.print_logger('\n', detector.log)

        # Long event search
        if not skglow:
            sm.print_logger('Starting search for glows...', detector.log)
            # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
            times = long_event_cutoff(detector)

            # Calling the long event search algorithm
            long_event_search(detector, times)

            sm.print_logger('Done.', detector.log)
            sm.print_logger('\n', detector.log)

    # Low memory mode:
    except MemoryError:
        sm.print_logger('\n', detector.log)
        sm.print_logger('Not enough memory. Entering low memory mode...', detector.log)
        sm.print_logger('\n', detector.log)
        # Measures the total combined size of all the data files
        total_file_size = 0
        for scint in detector:
            master_filelist = detector.attribute_retriever(scint, 'filelist')
            # Clears leftover data (just to be sure)
            detector.attribute_updator(scint, ['time', 'energy', 'filetime_extrema'],
                                       [np.array([]), np.array([]), []])
            for file in master_filelist:
                total_file_size += os.path.getsize(file)

        gc.collect()

        # Determines the appropriate number of chunks to split the day into
        operating_memory = sm.memory_allowance()
        available_memory = psutil.virtual_memory()[1]/4
        allowed_memory = available_memory

        num_chunks = 2  # minimum number of chunks to split the day into
        max_chunks = 16
        while num_chunks < max_chunks:
            mem_per_chunk = total_file_size/num_chunks
            if allowed_memory/(mem_per_chunk + operating_memory) >= 1:
                break

            num_chunks += 1

        # Makes the chunks
        chunk_list = []

        for chunk_num in range(1, num_chunks + 1):
            chunk = sc.Chunk(unit, first_sec, modes)
            chunk.log = log
            chunk_list.append(chunk)

        chunk_scint_list = chunk_list[0].scint_list

        for scint in detector:
            current_filelist = detector.attribute_retriever(scint, 'filelist')
            filelist_len = len(current_filelist)
            chunk_num = 1
            for chunk in chunk_list:
                if chunk_num == num_chunks:
                    chunk.attribute_updator(scint, 'filelist', current_filelist)
                else:
                    filelist_chunk = current_filelist[:int(filelist_len/num_chunks)]
                    chunk.attribute_updator(scint, 'filelist', filelist_chunk)
                    current_filelist = current_filelist[int(filelist_len/num_chunks):]

                chunk_num += 1

        # Imports data to each chunk and then pickles the chunks (and checks that data is actually present)
        print('Importing data...')
        print('\n')
        chunk_num = 1
        chunk_path_list = []
        # Temporary pickle feature for low memory mode. REMOVE WHEN PROGRAM IS FINISHED
        pickled_chunk_paths = glob.glob(f'{detector.results_loc}Results/{unit}/{full_day_string}/chunk*.pickle')
        pickled_chunk_paths.sort()
        if picklem and len(pickled_chunk_paths) > 0:
            chunk_path_list = pickled_chunk_paths
            data_present = True
        else:
            # Keeps timings consistent between chunks
            passtime = chunk_list[0].attribute_retriever(chunk_scint_list[0], 'passtime')
            passtime_dict = {scint: passtime.copy() for scint in chunk_scint_list}

            data_present = True
            for chunk in chunk_list:
                # Updates chunk to include previous chunk's passtime
                chunk.update_passtime(passtime_dict)

                sm.print_logger(f'Chunk {chunk_num} (of {num_chunks}):', detector.log)
                chunk.data_importer(existing_filelists=True)

                # Checking that data is present in the necessary scintillators
                if 'LP' in chunk.long_event_scint_list:
                    necessary_scintillators = chunk.long_event_scint_list
                else:
                    necessary_scintillators = chunk.long_event_scint_list + ['LP']

                for scint in necessary_scintillators:
                    if len(chunk.attribute_retriever(scint, 'time')) == 0:
                        data_present = False
                        print('\n\n')
                        print('\n', file=detector.log)
                        sm.print_logger('No/Missing data for specified day.', detector.log)
                        print('\n')
                        break

                if data_present:
                    print('\n\n')
                    # Makes a full list of filetime extrema for long event search
                    for scint in chunk:
                        extrema = detector.attribute_retriever(scint, 'filetime_extrema')
                        extrema += chunk.attribute_retriever(scint, 'filetime_extrema')
                        detector.attribute_updator(scint, 'filetime_extrema', extrema)

                    # Updates passtime
                    passtime_dict = chunk.return_passtime()

                    chunk_pickle_path = f'{detector.results_loc}Results/{unit}/{full_day_string}/' \
                                        f'chunk{chunk_num}.pickle'
                    chunk_path_list.append(chunk_pickle_path)
                    chunk_pickle = open(chunk_pickle_path, 'wb')
                    chunk.log = None
                    chunk.regex = None
                    pickle.dump(chunk, chunk_pickle)
                    chunk_pickle.close()

                    # Eliminates the chunk from active memory
                    del chunk
                    gc.collect()
                    chunk_list[chunk_num - 1] = 0

                    chunk_num += 1

                else:
                    break

        # Aborts the program for the day if necessary scintillator data is missing in any of the chunks
        if not data_present:
            for chunk_path in chunk_path_list:
                os.remove(chunk_path)

            break

        print('Done.')

        # Calibrates each scintillator
        if not skcali:
            print('\n')
            sm.print_logger('Calibrating scintillators and generating energy spectra...', detector.log)
            existing_spectra = {scint: np.array([]) for scint in chunk_scint_list}
            for chunk_path in chunk_path_list:
                chunk = sm.chunk_unpickler(chunk_path)
                chunk_spectra = chunk.make_spectra_hist(existing_spectra)
                del chunk

            # Calling the calibration algorithm
            detector.spectra_maker(existing_spectra=existing_spectra)

            sm.print_logger('Done.', detector.log)

        # Short event search
        if not skshort:
            sm.print_logger('\n', detector.log)
            sm.print_logger('Starting search for short events...', detector.log)
            sm.print_logger('\n', detector.log)
            existing_event_numbers = {scint: 0 for scint in chunk_scint_list}

            chunk_num = 1
            for chunk_path in chunk_path_list:
                chunk = sm.chunk_unpickler(chunk_path)
                chunk.log = log
                sm.print_logger(f'Chunk {chunk_num}:', detector.log)
                print('\n')

                # Calling the short event search algorithm
                existing_event_numbers = short_event_search(chunk, existing_event_numbers, low_mem=True)

                print('\n\n', file=detector.log)

                del chunk
                gc.collect()
                chunk_num += 1

            sm.print_logger('Done.', detector.log)
            sm.print_logger('\n', detector.log)

        # Long event search
        if not skglow:
            sm.print_logger('Starting search for glows...', detector.log)
            le_hist = np.array([])
            for chunk_path in chunk_path_list:
                chunk = sm.chunk_unpickler(chunk_path)
                # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
                times = long_event_cutoff(detector, chunk)

                # Histograms the counts from each chunk and combines them with the main one
                chunk_hist = long_event_search(detector, times, low_mem=True)
                le_hist = chunk_hist if len(le_hist) == 0 else le_hist + chunk_hist

                del chunk
                gc.collect()

            # Calling the long event search algorithm
            long_event_search(detector, np.array([]), existing_hist=le_hist)

            sm.print_logger('Done.', detector.log)
            sm.print_logger('\n', detector.log)

        log.close()
        # Deletes chunk .pickle files
        # REMOVE CONDITIONAL STATEMENT WHEN PROGRAM IS DONE
        if not picklem:
            for chunk_path in chunk_path_list:
                os.remove(chunk_path)

    # Missing data for necessary scintillators
    except FileNotFoundError:
        pass

    finally:
        del detector
        log.close()
        gc.collect()
