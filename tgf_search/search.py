"""A script that searches for TGFs and glows."""
import sys
import datetime as dt
import glob as glob
import os as os
import psutil as psutil
import gc as gc
import heapq
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit

import tgf_search.tools as tl
import tgf_search.parameters as params
from tgf_search.detectors.godot import Godot
from tgf_search.detectors.thor import Thor
from tgf_search.detectors.santis import Santis
from tgf_search.detectors.croatia import Croatia
from tgf_search.events.shortevent import ShortEvent
from tgf_search.events.longevent import LongEvent


# Returns the correct detector object based on the parameters provided
def get_detector(unit, date_str, modes=None, print_feedback=False):
    if modes is None:
        modes = list()

    if unit.upper() == 'GODOT':
        return Godot(unit, date_str, modes, print_feedback)
    elif unit.upper() == 'SANTIS':
        return Santis(unit, date_str, modes, print_feedback)
    elif unit.upper() == 'CROATIA':
        return Croatia(unit, date_str, modes, print_feedback)
    elif 'THOR' in unit.upper():
        if len(unit) >= 5 and unit[4].isnumeric() and int(unit[4]) <= 6:  # only 6 of them right now
            return Thor(unit, date_str, modes, print_feedback)
        else:
            raise ValueError(f"'{unit}' is not a valid detector.")
    else:
        raise ValueError(f"'{unit}' is not a valid detector.")


# Checks whether a short event is valid by passing it through several filters
def is_good_short_event(detector, modes, stats, times, energies, start, length):
    low_channel_counts = 0
    high_channel_counts = 0

    priority_queue = []

    stop = start + length

    clumpiness = 0
    clump_counts = 0

    # All filters share a single loop to speed things up
    high_channel_start = params.LOW_CHANNEL_START + params.CHANNEL_RANGE_WIDTH + params.CHANNEL_SEPARATION
    for i in range(start, stop):
        # Counting for low/high energy ratio filter
        if length >= params.GOOD_LEN_THRESH and not detector.is_named('THOR'):
            if params.LOW_CHANNEL_START <= energies[i] <= (params.LOW_CHANNEL_START + params.CHANNEL_RANGE_WIDTH):
                low_channel_counts += 1

            if high_channel_start <= energies[i] <= (high_channel_start + params.CHANNEL_RANGE_WIDTH):
                high_channel_counts += 1

        # Adding to priority queue for counts above minimum energy threshold filter
        heapq.heappush(priority_queue, -1 * energies[i])  # -1 to turn this into a max heap

        # Measuring clumpiness for successive crs filter
        if modes['aircraft'] and length < params.GOOD_LEN_THRESH and i > start:
            difference = times[i] - times[i - 1]
            if difference < params.DIFFERENCE_THRESH:
                clump_counts += 1
                if i == stop - 1:
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

    # Checks that there are fewer counts in the higher energy channels than the low energy channels
    # High/low channel ratio is not checked for THOR or for events with < GOOD_LEN_THRESH counts
    if not detector.is_named('THOR') and length >= params.GOOD_LEN_THRESH:
        if low_channel_counts == 0 or high_channel_counts / low_channel_counts > params.CHANNEL_RATIO:
            stats['removed_channel_ratio'] += 1
            return False

    # Eliminates the low energy events, which are likely just noise
    # At least min_noise_counts must be above noise_cutoff_energy
    for i in range(params.MIN_NOISE_COUNTS):
        if -1 * heapq.heappop(priority_queue) < params.NOISE_CUTOFF_ENERGY:
            stats['removed_low_energy'] += 1
            return False

    # Eliminates events that are just successive cosmic ray showers
    if modes['aircraft'] and length < params.GOOD_LEN_THRESH:
        clumpiness /= length
        if clumpiness >= params.CLUMPINESS_THRESH:
            stats['removed_crs'] += 1
            return False

    return True


# Calculates subscores, then uses them to calculate a final score for each event and then rank all the events
def rank_events(detector, potential_events, times, energies):
    weather_cache = {}
    for event in potential_events:
        event.calculate_score(weather_cache, detector, times, energies)

    ranked_events = sorted(potential_events, key=lambda x: -x.total_score)  # negative so we get a descending order sort
    for i in range(len(ranked_events)):
        ranked_events[i].rank = i + 1


# Short event search algorithm
def short_event_search(detector, modes, prev_event_numbers=None, low_mem=False):
    if modes['aircraft']:
        rollgap = params.AIRCRAFT_ROLLGAP
    elif modes['combo']:
        rollgap = params.COMBO_ROLLGAP
    else:
        rollgap = params.NORMAL_ROLLGAP

    event_numbers = prev_event_numbers if prev_event_numbers is not None else {}
    for scintillator in detector:
        if scintillator != detector.default_scintillator and not modes['allscints']:
            continue

        # Combining data from all available scintillators
        if modes['combo']:
            tl.print_logger('Searching combined scintillator data...', detector.log)
            (times,
             energies,
             wallclock,
             count_scints) = tl.combine_data(detector)
            gc.collect()

        # Normal operating mode (one scintillator at a time)
        else:
            tl.print_logger(f'Searching eRC {detector.get_attribute(scintillator, "eRC")} '
                            f'({scintillator})...', detector.log)
            times = detector.get_attribute(scintillator, 'time')
            energies = detector.get_attribute(scintillator, 'energy')
            wallclock = detector.get_attribute(scintillator, 'wc')
            count_scints = None

        stats = {
            'total_potential_events': 0,
            'total_threshold_reached': 0,
            'removed_len': 0,
            'removed_channel_ratio': 0,
            'removed_low_energy': 0,
            'removed_crs': 0
        }

        # Checks for an event by looking for a certain number of counts (rollgap + 1) in a small timeframe
        potential_events = []
        event_start = 0
        event_length = 0
        for i in range(len(times)):
            rolled_index = i + rollgap if i + rollgap < len(times) else (i + rollgap) - len(times)
            interval = times[rolled_index] - times[i]
            if 0 < interval <= params.SHORT_EVENT_TIME_SPACING:
                # Records the beginning index of a potential event
                if event_length == 0:
                    event_start = i
                    event_length = 1 + rollgap  # 1 for first count, rollgap for the others
                    stats['total_potential_events'] += 1
                # Measures the length of a potential event
                else:
                    event_length += 1

                # Counts the total number of times that the detection threshold was reached
                stats['total_threshold_reached'] += 1

            # Records the rough length of a potential event
            if interval > params.SHORT_EVENT_TIME_SPACING and event_length > 0:
                # Keeps potential event if it's longer than the specified minimum number of counts
                if event_length >= params.SHORT_EVENT_MIN_COUNTS:
                    # Runs potential event through filters
                    if is_good_short_event(detector, modes, stats,
                                           times, energies, event_start, event_length):
                        potential_events.append(ShortEvent(event_start, event_length, scintillator))

                else:
                    stats['removed_len'] += 1

                event_start = 0
                event_length = 0

        tl.print_logger(f'\n{len(potential_events)} potential events recorded', detector.log)
        if stats['removed_len'] > 0:
            tl.print_logger(f'{stats["removed_len"]} events removed due to insufficient length', detector.log)

        if stats['removed_channel_ratio'] > 0:
            tl.print_logger(f'{stats["removed_channel_ratio"]} events removed due to noise '
                            f'(bad high-to-low channel ratio)', detector.log)

        if stats['removed_low_energy'] > 0:
            tl.print_logger(f'{stats["removed_low_energy"]} events removed due to noise '
                            f'(minimum energy threshold not reached)', detector.log)

        if stats['removed_crs'] > 0:
            tl.print_logger(f'{stats["removed_crs"]} events removed due to successive CRS', detector.log)

        tl.print_logger(f'Detection threshold reached {stats["total_threshold_reached"]} times', detector.log)
        print('\n', file=detector.log)

        if len(potential_events) > 0:
            print('\n')
            tl.print_logger('Generating scatter plots and event files...', detector.log)
            print('\n', file=detector.log)

            if prev_event_numbers is not None:
                plots_already_made = event_numbers[scintillator]
            else:
                plots_already_made = 0

            assert plots_already_made <= params.MAX_PLOTS_PER_SCINT
            if (plots_already_made + len(potential_events)) >= params.MAX_PLOTS_PER_SCINT:
                max_plots = params.MAX_PLOTS_PER_SCINT - plots_already_made
            else:
                max_plots = len(potential_events)

            # Updates each event object to include its accurate subscores, total score, and rank
            rank_events(detector, potential_events, times, energies)

            plots_made = 0
            filecount_switch = True
            print('Potential short events:', file=detector.log)
            for i in range(len(potential_events)):
                # Stops making plots/event files/log entries after max reached
                # This is to prevent the program getting stuck on lab test days (with hundreds of thousands of "events")
                if (plots_made + plots_already_made) >= params.MAX_PLOTS_PER_SCINT:
                    print(f'Max number of loggable events ({params.MAX_PLOTS_PER_SCINT}) reached.',
                          file=detector.log)
                    break

                if not detector.gui:
                    print(f'{plots_made}/{max_plots}', end='\r')
                elif detector.gui and filecount_switch:
                    print(f'Making {max_plots} plots...')
                    filecount_switch = False

                event = potential_events[i]

                # Logs the event
                start_second = times[event.start] - params.SEC_PER_DAY if times[event.start] > params.SEC_PER_DAY else (
                    times)[event.start]
                print(f'{dt.datetime.utcfromtimestamp(times[event.start] + detector.first_sec)} UTC '
                      f'({start_second} seconds of day) - weather: {tl.weather_from_score(event.weather_subscore)}',
                      file=detector.log)
                print(f'    Rank: {event.rank}, Subscores: [Length: {event.len_subscore}, '
                      f'Clumpiness: {event.clumpiness_subscore}, HEL: {event.hel_subscore}]\n',
                      file=detector.log)

                if modes['combo']:
                    filelist_dict = dict()
                    filetime_extrema_dict = dict()
                    for j in range(event.start, event.stop):
                        if count_scints[j] not in filelist_dict:
                            filelist_dict[count_scints[j]] = detector.get_attribute(
                                count_scints[j], 'filelist')
                            filetime_extrema_dict[count_scints[j]] = detector.get_attribute(
                                count_scints[j], 'filetime_extrema')
                else:
                    filelist_dict = {scintillator: detector.get_attribute(scintillator, 'filelist')}
                    filetime_extrema_dict = {scintillator: detector.get_attribute(scintillator, 'filetime_extrema')}

                (event_file_dict,
                 filelist_dict,
                 filetime_extrema_dict) = event.get_filenames(filelist_dict, filetime_extrema_dict, times, count_scints)

                # Makes the scatter plot
                event.make_scatterplot(i + 1 + plots_already_made, event_file_dict, detector, times, energies,
                                       count_scints)

                # Makes the event file
                event.make_json(i + 1 + plots_already_made, event_file_dict, detector, times, energies,
                                wallclock, count_scints)

                plots_made += 1

            if not detector.gui:
                print(f'{plots_made}/{max_plots}\n', end='\r')

            event_numbers.update({scintillator: plots_made + plots_already_made})
        else:
            print('\n')

        # In combo mode, we only need to run through this loop once. The number of events found is tracked
        # in event_numbers by whichever scintillator the loop makes it to first (LP normally, NaI in 'allscints' mode)
        if modes['combo']:
            break

    if low_mem:
        return event_numbers


# Combines time data from several scintillators (if applicable) and cuts out counts below a certain energy
def long_event_cutoff(detector, modes, chunk=None):
    # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
    operating_obj = chunk if chunk is not None else detector  # operating_obj should contain data
    times = np.array([])
    energies = np.array([])
    # Checks to see that data is present in all preferred scintillators. Otherwise, uses the default scintillator.
    long_event_scintillators = operating_obj.long_event_scint_list
    for scintillator in operating_obj.long_event_scint_list:
        if not operating_obj.is_data_present(scintillator):
            long_event_scintillators = [operating_obj.default_scintillator]
            break

    for scintillator in long_event_scintillators:
        existing_calibration = True if len(
            detector.get_attribute(scintillator, 'calibration')) == 2 else False
        if scintillator == operating_obj.long_event_scint_list[0]:
            times = operating_obj.get_attribute(scintillator, 'time')
            if modes['skcali'] or not existing_calibration:
                energies = operating_obj.get_attribute(scintillator, 'energy')
            else:
                energies = tl.channel_to_mev(operating_obj.get_attribute(scintillator, 'energy'),
                                             detector.get_attribute(scintillator, 'calibration'),
                                             scintillator)
        else:
            times = np.append(times, operating_obj.get_attribute(scintillator, 'time'))
            if modes['skcali'] or not existing_calibration:
                energies = np.append(energies, operating_obj.get_attribute(scintillator, 'energy'))
            else:
                energies = np.append(energies,
                                     tl.channel_to_mev(operating_obj.get_attribute(scintillator, 'energy'),
                                                       detector.get_attribute(scintillator, 'calibration'),
                                                       scintillator))

    # Removes entries that are below a certain cutoff energy
    if (not detector.lp_calibrated and 'LP' in long_event_scintillators) or modes['skcali']:
        print('Missing calibration(s), beware radon washout!',
              file=detector.log)
    else:
        times = np.delete(times, np.where(energies < params.ENERGY_CUTOFF))

    if detector.processed:
        times = times - detector.first_sec

    return times


# Calculates mean for the long event search algorithm
def calculate_mue(hist_allday):
    hist_sum = 0
    hist_allday_nz = []
    for bin_val in hist_allday:
        if bin_val > 0:
            hist_allday_nz.append(bin_val)
            hist_sum += bin_val

    mue_val = hist_sum / len(hist_allday_nz)  # Mean number of counts/bin

    # Removes outliers (3 sigma) and recalculates mue
    abs_zscores = [np.abs((s - mue_val) / np.sqrt(mue_val)) for s in hist_allday_nz]
    sorting_order = np.argsort(abs_zscores)
    # Finds the index of the first zscore above 3 (this is essentially just a modified binary search algorithm)
    low = 0
    high = len(sorting_order) - 1
    while low < high:
        mid = low + (high - low) // 2
        if abs_zscores[sorting_order[mid]] >= 3:
            high = mid
        else:
            low = mid + 1

    hist_sum = 0
    for i in range(0, low):
        hist_sum += hist_allday_nz[sorting_order[i]]

    return hist_sum / low


# Calculates rolling mue/sigma baseline for the long event search algorithm
def calculate_rolling_baseline(day_bins, hist_allday, mue, sigma):
    center_index = 0
    num_bins = len(day_bins) - 1  # Hist array is shorter than bins array by 1 (-_-)
    so_poly = lambda t, a, b, c: a * t ** 2 + b * t + c
    line = lambda t, a, b: a * t + b

    l_bins = []
    l_counts = []
    l_bool = []
    l_zeros = 0

    r_bins = []
    r_counts = []
    r_bool = []
    r_zeros = 0

    # Setting up the initial right window
    for i in range(params.WINDOW_SIZE - 1):
        too_short = True if len(r_bool) == 0 else False
        index = center_index + params.WINDOW_GAP + 1 + i
        if hist_allday[index] > 0 and (too_short or hist_allday[i - 1]):
            r_bins.append(day_bins[index])
            r_counts.append(hist_allday[index])
            r_bool.append(True)
        else:
            r_zeros += 1
            r_bool.append(False)

    # Traversing the bins:
    while center_index < num_bins:
        # Advancing the left window
        # Adding to the front of the window
        l_index = center_index - (params.WINDOW_GAP + 1)
        if l_index >= 0:
            on_end = True if l_index - 1 < 0 else False
            # The bin being examined must be nonzero itself and have nonzero neighbors in order to be added
            # Neighbors are not checked for those bins on the ends, though
            if hist_allday[l_index] > 0 and (on_end or (
                    hist_allday[l_index - 1] > 0 and hist_allday[l_index + 1] > 0)):
                l_bins.append(day_bins[l_index])
                l_counts.append(hist_allday[l_index])
                l_bool.append(True)
            else:
                l_bool.append(False)
                l_zeros += 1

        # Removing from the back of the window
        if len(l_bool) > params.WINDOW_SIZE:
            if l_bool[0]:
                l_bins.pop(0)
                l_counts.pop(0)
            else:
                l_zeros -= 1

            l_bool.pop(0)

        # Advancing the right window
        # Adding to the front of the window
        r_index = center_index + params.WINDOW_GAP + params.WINDOW_SIZE
        if r_index < num_bins:
            on_end = True if r_index + 1 >= num_bins else False
            if hist_allday[r_index] > 0 and (on_end or (
                    hist_allday[r_index - 1] > 0 and hist_allday[r_index + 1] > 0)):
                r_bins.append(day_bins[r_index])
                r_counts.append(hist_allday[r_index])
                r_bool.append(True)
            else:
                r_bool.append(False)
                r_zeros += 1

        # Removing from the back of the window
        if len(r_bool) > 0 and center_index > 0:
            if r_bool[0]:
                r_bins.pop(0)
                r_counts.pop(0)
            else:
                r_zeros -= 1

            r_bool.pop(0)

        # Fitting curve is linear if either of the windows contains a zero value bin
        fit_curve = line if l_zeros != 0 or r_zeros != 0 else so_poly

        if len(l_bins) + len(r_bins) >= 3:
            fit, pcov = curve_fit(fit_curve, l_bins + r_bins, l_counts + r_counts)
            if fit_curve == so_poly:
                # Float casting is necessary for later parts of the day due to integer overflow
                mue[center_index] = so_poly(float(day_bins[center_index]), fit[0], fit[1], fit[2])
            else:
                mue[center_index] = line(day_bins[center_index], fit[0], fit[1])

            sigma[center_index] = np.sqrt(mue[center_index])

        center_index += 1

    return mue, sigma


# Long event search algorithm
def long_event_search(detector, modes, times, existing_hist=None, low_mem=False):
    # Makes one bin for every binsize seconds of the day (plus around 300 seconds more for the next day)
    day_bins = np.arange(0, 86700 + params.BIN_SIZE, params.BIN_SIZE)

    # Creates numerical values for histograms using numpy
    if existing_hist is not None:
        hist_allday = existing_hist
    else:
        hist_allday, bins_allday = np.histogram(times, bins=day_bins)
        if low_mem:
            return hist_allday

    # Calculates mean
    mue_val = calculate_mue(hist_allday)

    # Note: mue is an array full of the mue values for each bin. On a normal day, mue will be
    # filled entirely with mue_val. In aircraft mode, though, mue will be filled with a variety of different
    # values dictated by the rolling baseline algorithm in calculate_rolling_baseline
    mue = np.full(len(day_bins)-1, mue_val)
    sigma = np.full(len(day_bins)-1, np.sqrt(mue_val))

    # Calculating rolling baseline in aircraft mode
    if modes['aircraft']:
        mue, sigma = calculate_rolling_baseline(day_bins, hist_allday, mue, sigma)

    z_scores = []  # The z-scores themselves
    z_flags = []
    # Flags only those z-scores > params.FLAG_THRESH
    for i in range(len(hist_allday)):
        z_score = (hist_allday[i] - mue[i]) / sigma[i]
        z_scores.append(z_score)
        if z_score >= params.FLAG_THRESH:
            z_flags.append(i)

    # Sorts z-flags into actual potential glows
    potential_glows = []
    glow_start = 0
    glow_length = 0
    previous_time = 0
    for i in range(len(z_flags)):
        flag = z_flags[i]
        if glow_length == 0:  # First zscore only
            glow_start = flag
            glow_length += 1
        elif glow_length > 0 and day_bins[flag] - params.BIN_SIZE == previous_time:
            glow_length += 1
        elif (glow_length > 0 and day_bins[flag] - params.BIN_SIZE > previous_time) or i == len(z_flags) - 1:
            # Makes glow object
            potential_glows.append(LongEvent(glow_start, glow_length, z_scores, day_bins))
            glow_start = flag
            glow_length = 1

        previous_time = day_bins[flag]

    # Rejects events whose bins don't have enough counts (aircraft only)
    if modes['aircraft']:
        a_potential_glows = []
        for glow in potential_glows:
            for i in range(glow.start, glow.stop):
                if hist_allday[i] >= params.LONG_EVENT_MIN_COUNTS:
                    a_potential_glows.append(glow)
                    break

        potential_glows = a_potential_glows

    # Making histogram
    figure = plt.figure(figsize=[20, 11.0])
    plt.title(f'{detector.unit} at {detector.location["Location"]}, {detector.full_date_str}', loc='center')
    plt.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    ax1 = figure.add_subplot(5, 1, 1)  # Daily histogram
    # The 4 top events in descending order (if they exist)
    ax2 = figure.add_subplot(5, 1, 2)
    ax3 = figure.add_subplot(5, 1, 3)
    ax4 = figure.add_subplot(5, 1, 4)
    ax5 = figure.add_subplot(5, 1, 5)

    ax1.bar(day_bins[:-1], hist_allday, alpha=0.5, color='r', width=params.BIN_SIZE)
    ax1.set_xlabel('Seconds of Day (UT)')
    ax1.set_ylabel('Counts/bin')
    ax1.plot(day_bins[:-1], mue + params.FLAG_THRESH * sigma, color='blue', linestyle='dashed')

    # Creates legend
    allday_data = mpatches.Patch(color='r', label='All Energies')
    allday_thresh_sigma = mpatches.Patch(color='blue', label=f'{params.FLAG_THRESH} Sigma Above All Energies',
                                         linestyle='dashed')
    ax1.legend(handles=[allday_data, allday_thresh_sigma, ], bbox_to_anchor=(1.05, 1), loc=1,
               borderaxespad=0.)
    ax1.grid(True)

    tl.print_logger('Done.', detector.log)

    # Making event files and subplots (if they exist)
    if len(potential_glows) == 0:
        tl.print_logger('\n', detector.log)
        tl.print_logger(f'There were no potential glows on {detector.full_date_str}', detector.log)
    else:
        # Logs potential glows and sorts them in descending order depending on their highest z-score
        tl.print_logger('\n', detector.log)
        tl.print_logger('Generating event files...', detector.log)
        highest_scores = []
        print('\n', file=detector.log)
        print('Potential glows:', file=detector.log)

        eventpath = f'{detector.get_results_loc()}/Results/{detector.unit}/{detector.date_str}/' \
                    f'event files/long events/'
        tl.make_path(eventpath)
        event_number = 1
        files_made = 0
        filecount_switch = True
        for i in range(len(potential_glows)):
            if not detector.gui:
                print(f'{files_made}/{len(potential_glows)}', end='\r')
            elif detector.gui and filecount_switch:
                print(f'Making {len(potential_glows)} event files...')
                filecount_switch = False

            glow = potential_glows[i]
            highest_score = glow.highest_score
            highest_scores.append(highest_score)
            print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + detector.first_sec)} UTC ({glow.start_sec} '
                  f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                  f'{highest_score}', file=detector.log)

            event_file = open(f'{eventpath}{detector.date_str}_event{event_number}_zscore'
                              f'{int(highest_score)}.txt', 'w')
            print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + detector.first_sec)} UTC ({glow.start_sec} '
                  f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                  f'{highest_score}', file=event_file)

            long_event_scintillators = detector.long_event_scint_list
            for scintillator in detector.long_event_scint_list:
                if not detector.is_data_present(scintillator):
                    long_event_scintillators = [detector.default_scintillator]
                    break

            for scintillator in long_event_scintillators:
                print(f'{scintillator}:', file=event_file)
                filelist = detector.get_attribute(scintillator, 'filelist')
                filetime_extrema = detector.get_attribute(scintillator, 'filetime_extrema')
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

        if not detector.gui:
            print(f'{files_made}/{len(potential_glows)}\n', end='\r')

        tl.print_logger('Done.', detector.log)

        potential_glows = [potential_glows[s] for s in np.argsort(highest_scores)[::-1]]

        # Makes the histogram subplots
        ax_list = [ax2, ax3, ax4, ax5]
        for i in range(4):
            try:
                glow = potential_glows[i]
                glow.make_hist_subplot(ax_list[i], day_bins, hist_allday, mue, sigma)
            except IndexError:
                continue

    plt.tight_layout()

    # Saves the histogram(s):
    tl.print_logger('\n', detector.log)
    tl.print_logger('Saving Histogram...', detector.log)
    hist_path = f'{detector.get_results_loc()}/Results/{detector.unit}/{detector.date_str}/'
    tl.make_path(hist_path)
    plt.savefig(f'{hist_path}{detector.date_str}_histogram.png', dpi=500)
    plt.close(figure)


def main():
    try:
        first_date = str(sys.argv[1])
        second_date = str(sys.argv[2])
        unit = str(sys.argv[3])
    except IndexError:
        print('Please provide a first date, a second date, and a unit name.')
        exit()

    try:
        mode_info = sys.argv[4:]
    except IndexError:
        mode_info = []

    modes = dict()
    # Modes for skipping over certain algorithms (mostly to speed up testing)
    modes['skcali'] = True if 'skcali' in mode_info else False  # Skip detector calibration
    modes['skshort'] = True if 'skshort' in mode_info else False  # Skip short event search
    modes['skglow'] = True if 'skglow' in mode_info else False  # SKip long event search

    # Aircraft mode
    modes['aircraft'] = True if 'aircraft' in mode_info else False

    # Pickle mode
    modes['pickle'] = True if 'pickle' in mode_info else False

    # All scintillators mode (all the scintillators will be checked by the short event search algorithm)
    modes['allscints'] = True if 'allscints' in mode_info else False

    # Combo mode (all scintillator data is combined into one set of arrays and examined by the short event search algo)
    modes['combo'] = True if 'combo' in mode_info else False

    # Makes a list of all the dates on the requested range
    if int(second_date) < int(first_date):
        print('Not a valid date range.')
        exit()

    requested_dates = tl.make_date_list(first_date, second_date)
    for date_str in requested_dates:
        low_memory_mode = False
        day = int(date_str[4:])
        month = int(date_str[2:4])
        year = int('20' + date_str[0:2])
        full_day_str = dt.date(year, month, day)  # In format yyyy-mm-dd
        print(f'\n{full_day_str}:')

        # Initializes the detector object
        print('Importing data...')
        try:
            detector = get_detector(unit, date_str, mode_info, print_feedback=True)
        except ValueError:
            print('Not a valid detector.')
            exit()

        # Logs relevant data files and events in a .txt File
        log_path = f'{detector.get_results_loc()}/Results/{unit}/{date_str}/'
        tl.make_path(log_path)
        log = open(f'{log_path}log.txt', 'w')

        # Normal operating mode
        try:
            # Imports the data
            if modes['pickle']:  # In pickle mode: reads the pickle file, or exports one if it doesn't exist yet
                path_form = f'{detector.get_results_loc()}/Results/{detector.unit}/{detector.date_str}/detector.pickle'
                pickle_paths = glob.glob(path_form)
                if len(pickle_paths) > 0:
                    detector = tl.unpickle_detector(pickle_paths[0], mode_info)
                    detector.log = log
                else:
                    detector.log = log
                    detector.import_data()
                    tl.pickle_detector(detector, path_form)
                    detector.log = log
            else:
                detector.log = log
                detector.import_data()

            # raise MemoryError  # for low memory mode testing

            # Checks to see if there is actually data for the day
            if not detector:
                print('\n\n')
                print('\n', file=detector.log)
                tl.print_logger('No/Missing data for specified day.\n', detector.log)
                print('\n')
            # Otherwise runs normally
            else:
                print('\n\n')
                print('Done.')

                # Calibrates each scintillator
                if not modes['skcali']:
                    print('\n')
                    tl.print_logger('Calibrating scintillators and generating energy spectra...', detector.log)

                    # Calling the calibration algorithm
                    detector.calibrate(plot_spectra=True)

                    tl.print_logger('Done.', detector.log)

                # Short event search
                if not modes['skshort']:
                    tl.print_logger('\n', detector.log)
                    tl.print_logger('Starting search for short events...', detector.log)
                    tl.print_logger('\n', detector.log)

                    # Calling the short event search algorithm
                    short_event_search(detector, modes)

                    tl.print_logger('Done.', detector.log)
                    tl.print_logger('\n', detector.log)

                # Long event search
                if not modes['skglow']:
                    tl.print_logger('Starting search for glows...', detector.log)
                    # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
                    le_times = long_event_cutoff(detector, modes)

                    # Calling the long event search algorithm
                    long_event_search(detector, modes, le_times)

                    tl.print_logger('Done.', detector.log)
                    tl.print_logger('\n', detector.log)

        except MemoryError:
            low_memory_mode = True

        except FileNotFoundError:  # Missing necessary data
            pass

        # Low memory mode
        if low_memory_mode:
            try:
                tl.print_logger('\n', detector.log)
                tl.print_logger('Not enough memory. Entering low memory mode...', detector.log)
                tl.print_logger('\n', detector.log)
                # Measures the total combined size of all the data files
                total_file_size = 0
                for scintillator in detector:
                    master_filelist = detector.get_attribute(scintillator, 'filelist')
                    # Clears leftover data (just to be sure)
                    detector.set_attribute(scintillator, ['time', 'energy', 'filetime_extrema'],
                                              [np.array([]), np.array([]), []])
                    for file in master_filelist:
                        total_file_size += os.path.getsize(file)

                gc.collect()

                # Determines the appropriate number of chunks to split the day into
                allowed_memory = psutil.virtual_memory()[1] * params.TOTAL_MEMORY_ALLOWANCE_FRAC
                num_chunks = 2  # minimum number of chunks to split the day into
                max_not_exceeded = False
                while num_chunks < params.MAX_CHUNKS:
                    mem_per_chunk = total_file_size / num_chunks
                    if allowed_memory / (mem_per_chunk + params.OPERATING_MEMORY_ALLOWANCE) >= 1:
                        max_not_exceeded = True
                        break

                    num_chunks += 1

                if not max_not_exceeded:
                    raise MemoryError('very low available memory on system.')

                # Makes the chunks
                chunk_list = []

                for chunk_num in range(1, num_chunks + 1):
                    chunk = get_detector(unit, date_str, mode_info, print_feedback=True)
                    chunk.log = log
                    chunk_list.append(chunk)

                chunk_scint_list = chunk_list[0].scint_list

                for scintillator in detector:
                    current_filelist = detector.get_attribute(scintillator, 'filelist')
                    filelist_len = len(current_filelist)
                    chunk_num = 1
                    for chunk in chunk_list:
                        if chunk_num == num_chunks:
                            chunk.set_attribute(scintillator, 'filelist', current_filelist)
                        else:
                            filelist_chunk = current_filelist[:int(filelist_len / num_chunks)]
                            chunk.set_attribute(scintillator, 'filelist', filelist_chunk)
                            current_filelist = current_filelist[int(filelist_len / num_chunks):]

                        chunk_num += 1

                # Imports data to each chunk and then pickles the chunks (and checks that data is actually present)
                print('Importing data...')
                print('\n')
                chunk_num = 1
                chunk_path_list = []
                # Temporary pickle feature for low memory mode. REMOVE WHEN PROGRAM IS FINISHED
                pickled_chunk_paths = glob.glob(f'{detector.get_results_loc()}/Results/{unit}/{date_str}/chunk*.pickle')
                pickled_chunk_paths.sort()
                missing_data = False
                if modes['pickle'] and len(pickled_chunk_paths) > 0:
                    chunk_path_list = pickled_chunk_paths
                else:
                    # Keeps timings consistent between chunks
                    passtime = chunk_list[0].get_attribute(chunk_scint_list[0], 'passtime')
                    passtime_dict = {scint: passtime.copy() for scint in chunk_scint_list}

                    for chunk in chunk_list:
                        # Updates chunk to include previous chunk's passtime
                        chunk.update_passtime(passtime_dict)

                        tl.print_logger(f'Chunk {chunk_num} (of {num_chunks}):', detector.log)
                        chunk.import_data(existing_filelists=True)

                        # Checking that data is present in the necessary scintillators

                        # Aborts the program for the day if necessary scintillator data is missing in any of the chunks
                        if not chunk:
                            missing_data = True
                            print('\n\n')
                            print('\n', file=detector.log)
                            tl.print_logger('No/Missing data for specified day.', detector.log)
                            print('\n')
                            for chunk_path in chunk_path_list:
                                os.remove(chunk_path)

                            break
                        # Otherwise runs normally
                        else:
                            print('\n\n')
                            # Makes a full list of filetime extrema for long event search
                            for scintillator in chunk:
                                extrema = detector.get_attribute(scintillator, 'filetime_extrema')
                                extrema += chunk.get_attribute(scintillator, 'filetime_extrema')
                                detector.set_attribute(scintillator, 'filetime_extrema', extrema)

                            # Updates passtime
                            passtime_dict = chunk.return_passtime()

                            chunk_path = (f'{detector.get_results_loc()}/Results/{unit}/{date_str}/'
                                          f'chunk{chunk_num}.pickle')
                            chunk_path_list.append(chunk_path)
                            tl.pickle_chunk(chunk, chunk_path)

                            # Eliminates the chunk from active memory
                            del chunk
                            gc.collect()
                            chunk_list[chunk_num - 1] = 0

                            chunk_num += 1

                if not missing_data:
                    print('Done.')

                    # Calibrates each scintillator
                    if not modes['skcali']:
                        print('\n')
                        tl.print_logger('Calibrating scintillators and generating energy spectra...', detector.log)
                        existing_spectra = {scintillator: np.array([]) for scintillator in chunk_scint_list}
                        for chunk_path in chunk_path_list:
                            chunk = tl.unpickle_chunk(chunk_path)
                            existing_spectra = chunk.make_spectra_hist(existing_spectra)
                            del chunk

                        # Calling the calibration algorithm
                        detector.calibrate(existing_spectra=existing_spectra, plot_spectra=True)

                        tl.print_logger('Done.', detector.log)

                    # Short event search
                    if not modes['skshort']:
                        tl.print_logger('\n', detector.log)
                        tl.print_logger('Starting search for short events...', detector.log)
                        tl.print_logger('Warning! In low memory mode, short events '
                                        'will be ranked on a per-chunk basis.',
                                        detector.log)
                        tl.print_logger('\n', detector.log)
                        existing_event_numbers = {scintillator: 0 for scintillator in chunk_scint_list}

                        chunk_num = 1
                        for chunk_path in chunk_path_list:
                            chunk = tl.unpickle_chunk(chunk_path)
                            chunk.log = log
                            tl.print_logger(f'Chunk {chunk_num}:', detector.log)
                            print('\n')

                            # Calling the short event search algorithm
                            existing_event_numbers = short_event_search(chunk, modes,
                                                                        existing_event_numbers, low_mem=True)

                            del chunk
                            gc.collect()
                            chunk_num += 1

                            tl.print_logger('Done.', detector.log)
                            tl.print_logger('\n', detector.log)

                    # Long event search
                    if not modes['skglow']:
                        tl.print_logger('Starting search for glows...', detector.log)
                        le_hist = np.array([])
                        for chunk_path in chunk_path_list:
                            chunk = tl.unpickle_chunk(chunk_path)
                            # Converts energy channels to MeV using the locations of peaks/edges
                            # obtained during calibration
                            le_times = long_event_cutoff(detector, modes, chunk)

                            # Histograms the counts from each chunk and combines them with the main one
                            chunk_hist = long_event_search(detector, modes, le_times, low_mem=True)
                            le_hist = chunk_hist if len(le_hist) == 0 else le_hist + chunk_hist

                            del chunk
                            gc.collect()

                        # Calling the long event search algorithm
                        long_event_search(detector, modes, np.array([]), existing_hist=le_hist)

                        tl.print_logger('Done.', detector.log)
                        tl.print_logger('\n', detector.log)

                    log.close()
                    # Deletes chunk .pickle files
                    # REMOVE CONDITIONAL STATEMENT WHEN PROGRAM IS DONE
                    if not modes['pickle']:
                        for chunk_path in chunk_path_list:
                            os.remove(chunk_path)

            except FileNotFoundError:
                pass

        del detector
        log.close()
        gc.collect()


if __name__ == '__main__':
    main()
