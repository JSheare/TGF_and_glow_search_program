import sys
import warnings
import datetime as dt
import pickle as pickle
import glob as glob
import os as os
import psutil as psutil
import gc as gc
import heapq
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
import search_classes as sc
import search_module as sm


# Combines time data from several scintillators (if applicable) and cuts out counts below a certain energy
def long_event_cutoff(detector_obj, chunk_obj=None):
    # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
    operating_obj = chunk_obj if chunk_obj is not None else detector_obj
    le_times = np.array([])
    le_energies = np.array([])
    # Checks to see that data is present in all preferred scintillators. Otherwise, defaults to the large plastic.
    long_event_scintillators = operating_obj.long_event_scint_list
    for scintillator in operating_obj.long_event_scint_list:
        if not operating_obj.is_data_present(scintillator):
            long_event_scintillators = [operating_obj.default_scintillator]
            break

    for scintillator in long_event_scintillators:
        existing_calibration = True if len(
            detector_obj.attribute_retriever(scintillator, 'calibration')) == 2 else False
        if scintillator == operating_obj.long_event_scint_list[0]:
            le_times = operating_obj.attribute_retriever(scintillator, 'time')
            if skcali or not existing_calibration:
                le_energies = operating_obj.attribute_retriever(scintillator, 'energy')
            else:
                le_energies = sm.channel_to_mev(operating_obj.attribute_retriever(scintillator, 'energy'),
                                                detector_obj.attribute_retriever(scintillator, 'calibration'),
                                                scintillator)
        else:
            le_times = np.append(le_times, operating_obj.attribute_retriever(scintillator, 'time'))
            if skcali or not existing_calibration:
                le_energies = np.append(le_energies, operating_obj.attribute_retriever(scintillator, 'energy'))
            else:
                le_energies = np.append(le_energies,
                                        sm.channel_to_mev(operating_obj.attribute_retriever(scintillator, 'energy'),
                                                          detector_obj.attribute_retriever(scintillator, 'calibration'),
                                                          scintillator))

    # Removes entries that are below a certain cutoff energy
    if (not detector_obj.good_lp_calibration and 'LP' in long_event_scintillators) or skcali:
        print('Missing calibration(s), beware radon washout!',
              file=detector_obj.log)
    else:
        energy_cutoff = 1.9  # MeV
        le_times = np.delete(le_times, np.where(le_energies < energy_cutoff))

    if detector_obj.processed:
        le_times = le_times - first_sec

    return le_times


# Checks whether a short event is valid by passing it through several filters
def is_good_short_event(detector_obj, stats, times, energies, start, length):
    # Noise filter parameters:
    channel_range_width = 300
    channel_ratio = 0.5
    min_noise_counts = 3
    noise_cutoff_energy = 300

    # Successive crs filter parameters:
    difference_threshold = 2e-6
    gap_threshold = 10e-6
    clumpiness_threshold = 0.27

    low_channel_counts = 0
    high_channel_counts = 0

    priority_queue = []

    stop = start + length

    clumpiness = 0
    clump_counts = 0

    # All filters share a single loop to speed things up
    for i in range(start, stop):
        # Counting for low/high energy ratio filter
        if length >= 30 and not detector_obj.THOR:
            if 200 <= energies[i] <= (200 + channel_range_width):
                low_channel_counts += 1

            if (300 + channel_range_width) <= energies[i] <= (300 + 2 * channel_range_width):
                high_channel_counts += 1

        # Adding to priority queue for counts above minimum energy threshold filter
        heapq.heappush(priority_queue, -1 * energies[i])  # -1 to turn this into a max heap

        # Measuring clumpiness for successive crs filter
        if aircraft and length < 30 and i > start:
            difference = times[i] - times[i - 1]
            if difference < difference_threshold:
                clump_counts += 1
                if i == stop - 1:
                    clump_counts += 1

            else:
                # Adding to clumpiness when there's a clump of three or more
                if clump_counts >= 3:
                    clumpiness += 1
                    clump_counts = 0

                # Adding to the clumpiness when the gap between sufficient counts is greater than the threshold
                if difference >= gap_threshold and clump_counts == 0:
                    clumpiness += 1

                clump_counts = 0

    # Checks that there are fewer counts in the higher energy channels than the low energy channels
    # High/low channel ratio is not checked for THOR or for events with < 30 counts
    if not detector_obj.THOR and length >= 30:
        if low_channel_counts == 0 or high_channel_counts / low_channel_counts > channel_ratio:
            stats['removed_channel_ratio'] += 1
            return False

    # Eliminates the low energy events, which are likely just noise
    # At least min_noise_counts must be above noise_cutoff_energy
    for i in range(min_noise_counts):
        if -1 * heapq.heappop(priority_queue) < noise_cutoff_energy:
            stats['removed_low_energy'] += 1
            return False

    # Eliminates events that are just successive cosmic ray showers
    if aircraft and length < 30:
        clumpiness /= length
        if clumpiness >= clumpiness_threshold:
            stats['removed_crs'] += 1
            return False

    return True


# Calculates the ranking subscores for a potential short event
def calculate_subscores(detector_obj, event, times, energies, weather_cache):
    # Length parameters
    max_score_len = 30

    # Clumpiness parameters
    difference_threshold = 2e-6
    gap_threshold = 10e-6
    clumpiness_threshold = 0.27
    tossup = 0.2

    # High energy lead parameters
    energy_thresh = 15000

    clumpiness = 0
    clump_counts = 0

    high_energy_lead = 0
    leading_counts = 0

    for i in range(event.start, event.stop):
        if i > event.start:
            difference = times[i] - times[i - 1]
            if difference < difference_threshold:
                clump_counts += 1
                if clump_counts == 1:
                    leading_counts += 1
                    if energies[i - 1] >= energy_thresh:  # Clump starts on the previous index
                        high_energy_lead += 1

                # For the last count in the event
                if i == event.stop - 1:
                    clump_counts += 1

            else:
                # Adding to clumpiness when there's a clump of three or more
                if clump_counts >= 3:
                    clumpiness += 1
                    clump_counts = 0

                # Adding to the clumpiness when the gap between sufficient counts is greater than the threshold
                if difference >= gap_threshold and clump_counts == 0:
                    clumpiness += 1

                clump_counts = 0

    clumpiness /= event.length
    if leading_counts > 0:
        high_energy_lead /= leading_counts

    # Calculating the length subscore
    len_subscore = 1/max_score_len * event.length

    # Calculating the clumpiness subscore
    clumpiness_subscore = 1/(1 + np.e**(40*(clumpiness - tossup)))

    # Calculaing high energy leading count subscore
    hel_subscore = np.e ** -(5.3 * high_energy_lead)

    # Getting weather subscore
    if detector_obj.location['Nearest weather station'] != '':
        local_date, local_time = sm.convert_to_local(detector_obj.full_date_str, times[event.start])
        weather_subscore = sm.get_weather_conditions(local_date, local_time, detector_obj, weather_cache)
    else:
        weather_subscore = -1

    return len_subscore, clumpiness_subscore, hel_subscore, weather_subscore


# Calculates subscores, then uses them to calculate a final score for each event and then rank all the events
def rank_events(detector_obj, potential_events, times, energies):
    # Subscore weights
    len_weight = 0.3
    clumpiness_weight = 0.2
    hel_weight = 0.2
    weather_weight = 0.3
    assert len_weight + clumpiness_weight + hel_weight + weather_weight <= 1

    weather_cache = {}
    subscores = []
    scores = []
    for event in potential_events:
        length, clumpiness, hel, weather = calculate_subscores(detector_obj, event, times, energies, weather_cache)
        subscores.append((length, clumpiness, hel, weather))
        # If weather info couldn't be obtained, the weather subscore is removed and the remaining weights are adjusted
        # so that they stay proportional to one another
        if weather == -1:
            proportionality = 1/(len_weight + clumpiness_weight + hel_weight)
            score = proportionality * (len_weight * length + clumpiness_weight * clumpiness + hel_weight * hel)
        else:
            score = len_weight * length + clumpiness_weight * clumpiness + hel_weight * hel + weather_weight * weather

        scores.append(-score)  # Negative so that we get a descending order sort

    rankings = np.argsort(scores) + 1
    return [(event, subscores, rank) for event, subscores, rank in zip(potential_events, subscores, rankings)]


# Short event search algorithm
def short_event_search(detector_obj, prev_event_numbers=None, low_mem=False):
    # Parameters:
    rollgap = 18 if aircraft else 4
    event_time_spacing = 1e-3  # 1 millisecond
    event_min_counts = 10

    event_numbers = prev_event_numbers if prev_event_numbers is not None else {}
    for scintillator in detector_obj:
        if scintillator != 'LP' and not allscints:
            continue

        sm.print_logger(f'Searching eRC {detector_obj.attribute_retriever(scintillator, "eRC")} ({scintillator})...',
                        detector_obj.log)
        times = detector_obj.attribute_retriever(scintillator, 'time')
        energies = detector_obj.attribute_retriever(scintillator, 'energy')
        wallclock = detector_obj.attribute_retriever(scintillator, 'wc')
        filelist = detector_obj.attribute_retriever(scintillator, 'filelist')
        filetime_extrema = detector_obj.attribute_retriever(scintillator, 'filetime_extrema')

        stats = {
            'total_potential_events':   0,
            'total_threshold_reached':  0,
            'removed_len':              0,
            'removed_channel_ratio':    0,
            'removed_low_energy':       0,
            'removed_crs':              0
        }

        # Checks for an event by looking for a certain number of counts (rollgap + 1) in a small timeframe
        potential_event_list = []
        event_start = 0
        event_length = 0
        for i in range(len(times)):
            rolled_index = i + rollgap if i + rollgap < len(times) else (i + rollgap) - len(times)
            interval = times[rolled_index] - times[i]
            if 0 < interval <= event_time_spacing:
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
            if interval > event_time_spacing and event_length > 0:
                # Keeps potential event if it's longer than the specified minimum number of counts
                if event_length >= event_min_counts:
                    # Runs potential event through filters
                    if is_good_short_event(detector_obj, stats, times, energies, event_start, event_length):
                        potential_event_list.append(sc.ShortEvent(event_start, event_length, scintillator))

                else:
                    stats['removed_len'] += 1

                event_start = 0
                event_length = 0

        sm.print_logger(f'\n{len(potential_event_list)} potential events recorded', detector_obj.log)
        if stats['removed_len'] > 0:
            sm.print_logger(f'{stats["removed_len"]} events removed due to insufficient length', detector_obj.log)

        if stats['removed_channel_ratio'] > 0:
            sm.print_logger(f'{stats["removed_channel_ratio"]} events removed due to noise '
                            f'(bad high-to-low channel ratio)', detector_obj.log)

        if stats['removed_low_energy'] > 0:
            sm.print_logger(f'{stats["removed_low_energy"]} events removed due to noise '
                            f'(minimum energy threshold not reached)', detector_obj.log)

        if stats['removed_crs'] > 0:
            sm.print_logger(f'{stats["removed_crs"]} events removed due to successive CRS', detector_obj.log)

        sm.print_logger(f'Detection threshold reached {stats["total_threshold_reached"]} times', detector_obj.log)
        print('\n', file=detector_obj.log)

        if len(potential_event_list) > 0:
            print('\n')
            sm.print_logger('Generating scatter plots and event files...', detector_obj.log)
            print('\n', file=detector_obj.log)

            # Subplot timescales
            ts1 = 1e-4   # 100 microseconds
            ts2 = 0.005  # 5 milliseconds
            ts3 = 2      # 100 milliseconds
            ts_list = [ts1, ts2, ts3]

            if prev_event_numbers is not None:
                plots_already_made = event_numbers[scintillator]
            else:
                plots_already_made = 0

            max_plots_total = 1000 if not aircraft else float('inf')
            assert plots_already_made <= max_plots_total
            if (plots_already_made + len(potential_event_list)) >= max_plots_total:
                max_plots = max_plots_total - plots_already_made
            else:
                max_plots = len(potential_event_list)

            # Remakes the event list to include subscores and ranking for each event
            potential_event_list = rank_events(detector_obj, potential_event_list, times, energies)

            plots_made = 0
            filecount_switch = True
            print('Potential short events:', file=detector_obj.log)
            for i in range(len(potential_event_list)):
                # Stops making plots/event files/log entries after max reached
                # This is to prevent the program getting stuck on lab test days (with hundreds of thousands of "events")
                if (plots_made + plots_already_made) >= max_plots_total:
                    print(f'Max number of loggable events ({max_plots_total}) reached.', file=detector_obj.log)
                    break

                if not detector_obj.GUI:
                    print(f'{plots_made}/{max_plots}', end='\r')
                elif detector_obj.GUI and filecount_switch:
                    print(f'Making {max_plots} plots...')
                    filecount_switch = False

                event, subscores, rank = potential_event_list[i]
                weather_score = subscores[3]

                # Logs the event
                start_second = times[event.start] - 86400 if times[event.start] > 86400 else times[event.start]
                print(f'{dt.datetime.utcfromtimestamp(times[event.start] + first_sec)} UTC '
                      f'({start_second} seconds of day) - weather: {sm.weather_from_score(weather_score)}',
                      file=detector_obj.log)
                print(f'    Rank: {rank}, '
                      f'Subscores: [Length: {subscores[0]}, Clumpiness: {subscores[1]}, HEL: {subscores[2]}]\n',
                      file=detector_obj.log)

                # Note: I know the parameter lists below are long, but I'd rather have some ugly syntax here
                # than waste memory calling attribute_retriever and making unnecessary function-scoped copies

                # Makes the scatter plot
                event_file, filelist, filetime_extrema = event.get_filename(times, filelist, filetime_extrema)
                event.scatterplot_maker(ts_list, detector_obj, times, energies,
                                        i + 1 + plots_already_made, event_file, weather_score, rank)

                # Makes the event file
                event.json_maker(detector_obj, times, energies, wallclock, i + 1 + plots_already_made, event_file, rank)

                plots_made += 1

            if not detector_obj.GUI:
                print(f'{plots_made}/{max_plots}\n', end='\r')

            event_numbers.update({scintillator: plots_made + plots_already_made})
        else:
            print('\n')

    if low_mem:
        return event_numbers


# Long event search algorithm
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
    hist_sum = 0
    hist_allday_nz = []
    for bin_val in hist_allday:
        if bin_val > 0:
            hist_allday_nz.append(bin_val)
            hist_sum += bin_val

    mue_val = hist_sum / len(hist_allday_nz)  # Mean number of counts/bin

    # Removes outliers (3 sigma) and recalculates mue
    abs_zscores = [np.abs((s - mue_val)/np.sqrt(mue_val)) for s in hist_allday_nz]
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

    mue_val = hist_sum / len(hist_allday_nz)

    mue = np.full(len(day_bins)-1, mue_val)
    sigma = np.full(len(day_bins)-1, np.sqrt(mue_val))

    if aircraft:
        center_index = 0
        window_size = 20
        gap = 5
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
        for i in range(window_size - 1):
            too_short = True if len(r_bool) == 0 else False
            index = center_index + gap + 1 + i
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
            l_index = center_index - (gap + 1)
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
            if len(l_bool) > window_size:
                if l_bool[0]:
                    l_bins.pop(0)
                    l_counts.pop(0)
                else:
                    l_zeros -= 1

                l_bool.pop(0)

            # Advancing the right window
            # Adding to the front of the window
            r_index = center_index + gap + window_size
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
                params, pcov = curve_fit(fit_curve, l_bins + r_bins, l_counts + r_counts)
                if fit_curve == so_poly:
                    # Float casting is necessary for later parts of the day due to integer overflow
                    mue[center_index] = so_poly(float(day_bins[center_index]), params[0], params[1], params[2])
                else:
                    mue[center_index] = line(day_bins[center_index], params[0], params[1])

                sigma[center_index] = np.sqrt(mue[center_index])

            center_index += 1

    z_scores = []  # The z-scores themselves
    z_flags = []
    # Flags only those z-scores > flag_threshold
    for i in range(len(hist_allday)):
        z_score = (hist_allday[i] - mue[i]) / sigma[i]
        z_scores.append(z_score)
        if z_score >= flag_threshold:
            z_flags.append(i)

    # Sorts z-flags into actual potential glows
    glow_start = 0
    glow_length = 0
    previous_time = 0
    potential_glow_list = []
    for i in range(len(z_flags)):
        flag = z_flags[i]
        if glow_length == 0:  # First zscore only
            glow_start = flag
            glow_length += 1
        elif glow_length > 0 and day_bins[flag] - binsize == previous_time:
            glow_length += 1
        elif (glow_length > 0 and day_bins[flag] - binsize > previous_time) or i == len(z_flags) - 1:
            # Makes glow object
            potential_glow_list.append(sc.PotentialGlow(glow_start, glow_length, z_scores, day_bins, binsize))
            glow_start = flag
            glow_length = 1

        previous_time = day_bins[flag]

    # Rejects events whose bins don't have enough counts (aircraft only)
    if aircraft:
        a_potential_glow_list = []
        minimum_counts = 1000
        for glow in potential_glow_list:
            for i in range(glow.start, glow.stop):
                if hist_allday[i] > minimum_counts:
                    a_potential_glow_list.append(glow)
                    break

        potential_glow_list = a_potential_glow_list

    sm.print_logger('Done.', detector_obj.log)
    if len(potential_glow_list) == 0:
        sm.print_logger('\n', detector_obj.log)
        sm.print_logger(f'There were no potential glows on {full_day_str}', detector_obj.log)
    else:
        # Logs potential glows and sorts them in descending order depending on their highest z-score
        sm.print_logger('\n', detector_obj.log)
        sm.print_logger('Generating event files...', detector_obj.log)
        highest_scores = []
        print('\n', file=detector_obj.log)
        print('Potential glows:', file=detector_obj.log)

        eventpath = f'{detector_obj.results_loc}Results/{detector_obj.unit}/{detector_obj.date_str}/' \
                    f'event files/long events/'
        sm.path_maker(eventpath)
        event_number = 1
        files_made = 0
        filecount_switch = True
        for i in range(len(potential_glow_list)):
            if not detector_obj.GUI:
                print(f'{files_made}/{len(potential_glow_list)}', end='\r')
            elif detector_obj.GUI and filecount_switch:
                print(f'Making {len(potential_glow_list)} event files...')
                filecount_switch = False

            glow = potential_glow_list[i]
            highest_score = glow.highest_score
            highest_scores.append(highest_score)
            print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + first_sec)} UTC ({glow.start_sec} '
                  f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                  f'{highest_score}', file=detector_obj.log)

            event_file = open(f'{eventpath}{detector_obj.date_str}_event{event_number}_zscore'
                              f'{int(highest_score)}.txt', 'w')
            print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + first_sec)} UTC ({glow.start_sec} '
                  f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                  f'{highest_score}', file=event_file)

            long_event_scintillators = detector_obj.long_event_scint_list
            for scintillator in detector_obj.long_event_scint_list:
                if not detector_obj.is_data_present(scintillator):
                    long_event_scintillators = [detector_obj.default_scintillator]
                    break

            for scintillator in long_event_scintillators:
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

        if not detector_obj.GUI:
            print(f'{files_made}/{len(potential_glow_list)}\n', end='\r')

        sm.print_logger('Done.', detector_obj.log)

        potential_glow_list = [potential_glow_list[s] for s in np.argsort(highest_scores)[::-1]]

        # Plotting the histograms
        sm.print_logger('\n', detector_obj.log)
        sm.print_logger('Generating Histogram...', detector_obj.log)
        figu = plt.figure(figsize=[20, 11.0])
        plt.title(f'{unit} at {detector_obj.location["Location"]}, {str(full_day_str)}', loc='center')
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

        # Saves the histograms:
        hist_path = f'{detector_obj.results_loc}Results/{unit}/{date_str}/'
        sm.path_maker(hist_path)
        plt.savefig(f'{hist_path}{date_str}_histogram.png', dpi=500)
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
    date_str = first_date
    century = '20'
    while True:
        date_str = sm.roll_date_forward(date_str)

        if date_str == second_date:
            requested_dates.append(date_str)
            break
        else:
            requested_dates.append(date_str)

for date in requested_dates:
    low_memory_mode = False
    date_str = str(date)  # In format yymmdd
    day = int(date_str[4:])
    month = int(date_str[2:4])
    year = int('20' + date_str[0:2])
    full_day_str = dt.date(year, month, day)  # In format yyyy-mm-dd
    print(f'\n{full_day_str}:')

    # EPOCH time conversions
    first_sec = (dt.datetime(year, month, day, 0, 0) - dt.datetime(1970, 1, 1)).total_seconds()
    print(f'The first second of {year}-{month}-{day} is: {int(first_sec)}')
    print(f'The last second of {year}-{month}-{day} is: {int(first_sec + 86400)}')

    # Initializes the detector object
    print('Importing data...')
    try:
        detector = sc.Detector(unit, first_sec, modes, print_feedback=True)
    except ValueError:
        print('Not a valid detector.')
        exit()

    # Logs relevant data files and events in a .txt File
    log_path = f'{detector.results_loc}Results/{unit}/{date_str}/'
    sm.path_maker(log_path)
    log = open(f'{log_path}log.txt', 'w')
    log.write(f'The first second of {year}-{month}-{day} is: {int(first_sec)}\n')
    log.write(f'The last second of {year}-{month}-{day} is: {int(first_sec + 86400)}\n')

    # Normal operating mode
    try:
        # Imports the data
        if picklem:
            pickle_path = glob.glob(f'{detector.results_loc}Results/{unit}/{date_str}/detector.pickle')
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
                detector_pickle = open(f'{detector.results_loc}Results/{unit}/{date_str}/detector.pickle',
                                       'wb')
                pickle.dump(detector, detector_pickle)
                detector_pickle.close()
                detector.log = log
        else:
            detector.log = log
            detector.data_importer()

        # Checks to see if there is actually data for the day
        if not detector:
            print('\n\n')
            print('\n', file=detector.log)
            sm.print_logger('No/Missing data for specified day.\n', detector.log)
            print('\n')
        # Otherwise runs normally
        else:
            print('\n\n')
            print('Done.')

            # Calibrates each scintillator
            if not skcali:
                print('\n')
                sm.print_logger('Calibrating scintillators and generating energy spectra...', detector.log)

                # Calling the calibration algorithm
                detector.calibrate()

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

    except MemoryError:
        low_memory_mode = True

    except FileNotFoundError:  # Missing necessary data
        pass

    # Low memory mode
    if low_memory_mode:
        try:
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
            available_memory = psutil.virtual_memory()[1] / 4
            allowed_memory = available_memory

            num_chunks = 2  # minimum number of chunks to split the day into
            max_chunks = 16
            max_not_exceeded = False
            while num_chunks < max_chunks:
                mem_per_chunk = total_file_size / num_chunks
                if allowed_memory / (mem_per_chunk + operating_memory) >= 1:
                    max_not_exceeded = True
                    break

                num_chunks += 1

            if not max_not_exceeded:
                raise MemoryError('MemoryError: very low available memory on system.')

            # Makes the chunks
            chunk_list = []

            for chunk_num in range(1, num_chunks + 1):
                chunk = sc.Chunk(unit, first_sec, modes, print_feedback=True)
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
                        filelist_chunk = current_filelist[:int(filelist_len / num_chunks)]
                        chunk.attribute_updator(scint, 'filelist', filelist_chunk)
                        current_filelist = current_filelist[int(filelist_len / num_chunks):]

                    chunk_num += 1

            # Imports data to each chunk and then pickles the chunks (and checks that data is actually present)
            print('Importing data...')
            print('\n')
            chunk_num = 1
            chunk_path_list = []
            # Temporary pickle feature for low memory mode. REMOVE WHEN PROGRAM IS FINISHED
            pickled_chunk_paths = glob.glob(f'{detector.results_loc}Results/{unit}/{date_str}/chunk*.pickle')
            pickled_chunk_paths.sort()
            missing_data = False
            if picklem and len(pickled_chunk_paths) > 0:
                chunk_path_list = pickled_chunk_paths
            else:
                # Keeps timings consistent between chunks
                passtime = chunk_list[0].attribute_retriever(chunk_scint_list[0], 'passtime')
                passtime_dict = {scint: passtime.copy() for scint in chunk_scint_list}

                for chunk in chunk_list:
                    # Updates chunk to include previous chunk's passtime
                    chunk.update_passtime(passtime_dict)

                    sm.print_logger(f'Chunk {chunk_num} (of {num_chunks}):', detector.log)
                    chunk.data_importer(existing_filelists=True)

                    # Checking that data is present in the necessary scintillators

                    # Aborts the program for the day if necessary scintillator data is missing in any of the chunks
                    if not chunk:
                        missing_data = True
                        print('\n\n')
                        print('\n', file=detector.log)
                        sm.print_logger('No/Missing data for specified day.', detector.log)
                        print('\n')
                        for chunk_path in chunk_path_list:
                            os.remove(chunk_path)

                        break
                    # Otherwise runs normally
                    else:
                        print('\n\n')
                        # Makes a full list of filetime extrema for long event search
                        for scint in chunk:
                            extrema = detector.attribute_retriever(scint, 'filetime_extrema')
                            extrema += chunk.attribute_retriever(scint, 'filetime_extrema')
                            detector.attribute_updator(scint, 'filetime_extrema', extrema)

                        # Updates passtime
                        passtime_dict = chunk.return_passtime()

                        chunk_pickle_path = f'{detector.results_loc}Results/{unit}/{date_str}/' \
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

            if not missing_data:
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
                    detector.calibrate(existing_spectra=existing_spectra)

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

        except FileNotFoundError:
            pass

    del detector
    log.close()
    gc.collect()
