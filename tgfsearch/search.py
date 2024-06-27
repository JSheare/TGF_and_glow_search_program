"""A script that searches for TGFs and glows."""
import sys as sys
import datetime as dt
import glob as glob
import os as os
import traceback as traceback
import psutil as psutil
import gc as gc
import heapq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit

# Adds parent directory to sys.path. Necessary to make the imports below work when running this file as a script
parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import tgfsearch.tools as tl
import tgfsearch.parameters as params
from tgfsearch.detectors.godot import Godot
from tgfsearch.detectors.thor import Thor
from tgfsearch.detectors.santis import Santis
from tgfsearch.detectors.croatia import Croatia
from tgfsearch.events.shortevent import ShortEvent
from tgfsearch.events.longevent import LongEvent


# A helper class used to keep track of aligned traces
class TraceInfo:
    def __init__(self, trace_name, buff_no, times, energies):
        self.trace_name = trace_name
        self.buff_no = buff_no
        self.times = times
        self.energies = energies


# Returns the correct detector object based on the parameters provided
def get_detector(unit, date_str, print_feedback=False):
    if unit.upper() == 'GODOT':
        return Godot(unit, date_str, print_feedback)
    elif unit.upper() == 'SANTIS':
        return Santis(unit, date_str, print_feedback)
    elif unit.upper() == 'CROATIA':
        return Croatia(unit, date_str, print_feedback)
    elif 'THOR' in unit.upper():
        if len(unit) >= 5 and unit[4].isnumeric() and int(unit[4]) <= 6:  # only 6 of them right now
            return Thor(unit, date_str, print_feedback)
        else:
            raise ValueError(f"'{unit}' is not a valid detector.")
    else:
        raise ValueError(f"'{unit}' is not a valid detector.")


# Makes the modes dict used by many of the program's functions
def get_modes(mode_info):
    modes = dict()
    # Modes for skipping over certain algorithms (mostly to speed up testing)
    modes['skcali'] = True if '--skcali' in mode_info else False  # Skip detector calibration
    modes['skshort'] = True if '--skshort' in mode_info else False  # Skip short event search
    modes['skglow'] = True if '--skglow' in mode_info else False  # SKip long event search

    # Aircraft mode
    modes['aircraft'] = True if '--aircraft' in mode_info else False

    # Pickle mode
    modes['pickle'] = True if '--pickle' in mode_info else False

    # All scintillators mode (all the scintillators will be checked by the short event search algorithm)
    modes['allscints'] = True if '--allscints' in mode_info else False

    # Combo mode (all scintillator data is combined into one set of arrays and examined by the short event search algo)
    modes['combo'] = True if '--combo' in mode_info else False

    # Template mode (make LP template)
    modes['template'] = True if '--template' in mode_info else False

    # GUI mode (running script from gui)
    modes['gui'] = True if '-g' in mode_info else False

    # Custom mode (use custom import and/or export directories
    modes['custom'] = True if '-c' in mode_info else False

    # Processed mode (use processed data). Only available for Godot
    modes['processed'] = True if '-p' in mode_info else False

    return modes


# Plots the given list of traces
def plot_traces(detector, scintillator, trace_names):
    if trace_names:
        # Makes the trace plot path
        plot_path = f'{detector.get_results_loc()}/traces'
        tl.make_path(plot_path)
        for trace_name in trace_names:
            trace = detector.get_trace(scintillator, trace_name, deepcopy=False)
            trace_single_buff = trace[trace['BufferNo'] == trace['BufferNo'].iloc[0]]
            plot_name = 'eRC' + trace_name.split('eRC')[-1][1:].split('.')[0]
            plt.plot(trace_single_buff['Seconds'], trace_single_buff['pulse'])
            plt.xlabel('Seconds')
            plt.ylabel('Pulse Magnitude')
            plt.title(plot_name)
            plt.savefig(f'{plot_path}/{plot_name}.png')
            plt.clf()


# Searches through all traces in detector and returns lists of good ones for each scintillator. Also plots good traces
def trace_search(detector):
    # Filtering traces and setting up trace_dict, which keeps track of filtered trace names
    trace_dict = {}
    for scintillator in detector:
        trace_dict[scintillator] = tl.filter_traces(detector, scintillator)

        # Making plots for the traces that passed through the filter
        plot_traces(detector, scintillator, trace_dict[scintillator])

    return trace_dict


# Checks whether a short event is valid by passing it through several filters
def is_good_short_event(detector, modes, stats, times, energies, start, length):
    # Checks that the length of the event is greater than or equal to a certain minimum number of counts
    if length < params.SHORT_EVENT_MIN_COUNTS:
        stats['removed_len'] += 1
        return False

    low_channel_counts = 0
    high_channel_counts = 0

    priority_queue = []

    stop = start + length

    clumpiness = 0
    clump_counts = 0

    # All filters share a single loop to speed things up
    high_channel_start = params.LOW_CHANNEL_START + params.CHANNEL_RANGE_WIDTH + params.CHANNEL_SEPARATION
    for i in range(start, stop):
        count_energy = energies[i]
        # Counting for low/high energy ratio filter
        if length >= params.GOOD_LEN_THRESH and not detector.is_named('THOR'):
            if params.LOW_CHANNEL_START <= count_energy <= (params.LOW_CHANNEL_START + params.CHANNEL_RANGE_WIDTH):
                low_channel_counts += 1

            if high_channel_start <= count_energy <= (high_channel_start + params.CHANNEL_RANGE_WIDTH):
                high_channel_counts += 1

        # Adding to priority queue for counts above minimum energy threshold filter
        heapq.heappush(priority_queue, -1 * count_energy)  # -1 to turn this into a max heap

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
        event.calculate_score(detector, weather_cache, times, energies)

    ranked_events = sorted(potential_events, key=lambda x: -x.total_score)  # negative so we get a descending order sort
    for i in range(len(ranked_events)):
        ranked_events[i].rank = i + 1


# Short event search algorithm
def find_short_events(detector, modes, scintillator, rollgap, times, energies):
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
        elif interval > params.SHORT_EVENT_TIME_SPACING and event_length > 0:
            # Runs potential event through filters
            if is_good_short_event(detector, modes, stats, times, energies, event_start, event_length):
                potential_events.append(ShortEvent(event_start, event_length, scintillator))

            event_start = 0
            event_length = 0

    tl.print_logger(f'{len(potential_events)} potential events recorded', detector.log)
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

    return potential_events


# Finds the list mode file(s) associated with a short event
def find_se_files(detector, event, times, count_scints=None):
    event_lm_files = {}
    for i in range(event.start, event.stop):
        if count_scints is not None:
            scintillator = count_scints[i]
        else:
            scintillator = event.scintillator

        if scintillator not in event_lm_files:
            event_time = times[i]
            event_lm_files[scintillator] = detector.find_lm_file(scintillator, event_time)

        # So that we don't loop through the whole event for no reason when not in combo mode
        if count_scints is None:
            break

    return event_lm_files


# Finds the traces associated with a short event
def find_se_traces(detector, event, trace_dict, times, count_scints=None):
    event_traces = {}
    for i in range(event.start, event.stop):
        if count_scints is not None:
            scintillator = count_scints[i]
        else:
            scintillator = event.scintillator

        if scintillator not in event_traces:
            event_time = times[i]
            event_lm_file_data = detector.get_lm_file(scintillator, detector.find_lm_file(scintillator, event_time))
            potential_matches = detector.find_matching_traces(scintillator, event_time,
                                                              trace_list=trace_dict[scintillator])
            trace_info = None
            best_diff = float('inf')
            for trace_name in potential_matches:
                trace = detector.get_trace(scintillator, trace_name)
                buff_no = trace['BufferNo'].iloc[0]
                # Aligning the trace and checking how close it is to the event
                try:  # In case trace isn't in the same file after all
                    trace_times, trace_energies = tl.align_trace(trace, event_lm_file_data, buff_no)
                    diff = abs(trace_times[0] - event_time)
                    # The best match is whichever trace is closest to the event
                    if diff < best_diff:
                        best_diff = diff
                        trace_info = TraceInfo(trace_name, buff_no, trace_times, trace_energies)

                except ValueError:
                    continue

            if (trace_info is not None and len(np.where((trace_info.times >= times[event.start]) &
                                                        (trace_info.times <= times[event.stop]))[0]) > 0):
                event_traces[scintillator] = trace_info

            # So that we don't loop through the whole event for no reason when not in combo mode
            if count_scints is None:
                break

    return event_traces


# Makes the scatter plot for a short event
def make_se_scatterplot(detector, event, times, energies, count_scints):
    # Subplot timescales
    timescales = [params.SE_TIMESCALE_ONE, params.SE_TIMESCALE_TWO, params.SE_TIMESCALE_THREE]

    # Dot colors. Note that if an instrument with more than just NaI, SP, MP, and LP is ever added, this
    # will result in key errors
    colors = {'NaI': params.NAI_COLOR, 'SP': params.SP_COLOR, 'MP': params.MP_COLOR, 'LP': params.LP_COLOR}

    # Truncated time and energy arrays to speed up scatter plot making
    fraction_of_day = 1 / 128
    spacer = int((len(times) * fraction_of_day) / 2)
    left_edge = 0 if event.start - spacer < 0 else event.start - spacer
    right_edge = (len(times) - 1) if event.stop + spacer > (len(times) - 1) else event.stop + spacer
    if count_scints is not None:
        times_dict = tl.separate_data(times, count_scints, left_edge, right_edge)
        energies_dict = tl.separate_data(energies, count_scints, left_edge, right_edge)
    else:
        times_dict = {event.scintillator: times[left_edge:right_edge]}
        energies_dict = {event.scintillator: energies[left_edge:right_edge]}

    figure1 = plt.figure(figsize=[20, 11.0], dpi=150.)
    figure1.suptitle(f'{event.scintillator} Event {str(event.number)}, '
                     f'{dt.datetime.utcfromtimestamp(times[event.start] + detector.first_sec)} UTC, '
                     f'{event.length} counts \n Weather: {tl.weather_from_score(event.weather_subscore)} \n'
                     f'Rank: {event.rank}', fontsize=20)
    ax1 = figure1.add_subplot(3, 1, 1)
    ax2 = figure1.add_subplot(3, 1, 2)
    ax3 = figure1.add_subplot(3, 1, 3)
    ax_list = [ax1, ax2, ax3]
    assert len(ax_list) == len(timescales)

    event_times = times[event.start:event.stop]
    best_time = event_times[np.argmin(np.abs(event_times - np.roll(event_times, 1)))]
    for i in range(len(ax_list)):
        ts = timescales[i]
        ax = ax_list[i]
        ax.set_xlim(xmin=best_time - (ts / 2), xmax=best_time + (ts / 2))

        dot_size = 5 if ts == params.SE_TIMESCALE_THREE else 3  # makes larger dots for top plot
        ax.set_yscale('log')
        ax.set_ylim([0.6, 1e5])
        for scintillator in times_dict:
            # 0.6 to avoid annoying divide by zero warnings
            ax.scatter(times_dict[scintillator], energies_dict[scintillator] + 0.6,
                       s=dot_size, zorder=1, alpha=params.DOT_ALPHA, label=scintillator, color=colors[scintillator])

            # Plotting traces on only the first subplot
            if ts == params.SE_TIMESCALE_ONE and scintillator in event.traces:
                ax.plot(event.traces[scintillator].times, event.traces[scintillator].energies + 0.6,
                        zorder=-1, alpha=params.DOT_ALPHA, color=colors[scintillator])

        ax.set_xlabel(f'Time (Seconds, {ts}s total)')
        ax.set_ylabel('Energy Channel')
        # Lines appear (100*percent)% to the left or right of event start/stop depending on subplot timescale
        percent = 0.001
        # Event start line is orange, end line is red
        ax.vlines([event_times[0] - percent * ts, event_times[-1] + percent * ts], 0, 1e5,
                  colors=['orange', 'r'], linewidth=1, zorder=-1, alpha=0.3)

    # Adds a legend to the plot if we're in combo mode
    if count_scints is not None:
        plt.legend(loc='lower right')

    # Adds the name of the relevant list mode file (and trace, if applicable) to the scatter plot
    if count_scints is None:
        plt.title(f'Obtained from {event.lm_files[event.scintillator]}'
                  + (f'\nand {event.traces[event.scintillator].trace_name}' if
                     event.scintillator in event.traces else ''),
                  fontsize=15, y=-0.4)

    # Saves the scatter plot
    # Note: with this code, if an event happens in that 200-300 seconds of the next day that are included in the
    # last file, the image will have the wrong date in its name (though the timestamp in the scatter plot title will
    # always be correct)
    scatter_path = f'{detector.get_results_loc()}/scatter_plots'
    tl.make_path(scatter_path)
    event_num_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(event.number)))
    rank_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(event.rank)))
    figure1.savefig(f'{scatter_path}/{detector.date_str}_{event.scintillator}_'
                    f'event{event_num_padding}{event.number}_rank{rank_padding}{event.rank}.png')
    plt.close(figure1)


# Makes the json file for a short event
def make_se_json(detector, event, times, energies, wallclock, count_scints):
    eventpath = f'{detector.get_results_loc()}/event_files/short_events/'
    tl.make_path(eventpath)
    event_frame = pd.DataFrame()
    event_frame['wc'] = wallclock[event.start:event.stop]
    event_frame['SecondsOfDay'] = times[event.start:event.stop]
    event_frame['energies'] = energies[event.start:event.stop]
    event_frame['count_scintillator'] = count_scints[event.start:event.stop] if count_scints is not None else (
            [event.scintillator] * event.length)
    event_frame['lm_file'] = [event.lm_files[scintillator] for scintillator in event_frame['count_scintillator']]
    event_frame['trace_file'] = [(event.traces[scintillator].trace_name if scintillator in event.traces else '') for
                                 scintillator in event_frame['count_scintillator']]

    # Saves the json file
    event_num_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(event.number)))
    rank_padding = '0' * (len(str(params.MAX_PLOTS_PER_SCINT)) - len(str(event.rank)))
    event_frame.to_json(f'{eventpath}/{detector.date_str}_{event.scintillator}_'
                        f'event{event_num_padding}{event.number}_rank{rank_padding}{event.rank}.json')


# Runs the short event search
def short_event_search(detector, modes, trace_dict, prev_event_numbers=None):
    if modes['aircraft']:
        rollgap = params.AIRCRAFT_ROLLGAP
    elif modes['combo']:
        rollgap = params.COMBO_ROLLGAP
    else:
        rollgap = params.NORMAL_ROLLGAP

    event_numbers = prev_event_numbers if prev_event_numbers is not None else {}
    for i in range(len(detector.scint_list)):
        if detector.scint_list[i] != detector.default_scintillator and not modes['allscints']:
            continue

        tl.print_logger('', detector.log)
        # Combining data from all available scintillators (combo mode)
        if modes['combo']:
            scintillator = 'CM'
            tl.print_logger('Searching combined scintillator data...', detector.log)
            (times,
             energies,
             wallclock,
             count_scints) = tl.combine_data(detector)

        # Normal operating mode (one scintillator at a time)
        else:
            scintillator = detector.scint_list[i]
            tl.print_logger(f'Searching eRC {detector.get_attribute(scintillator, "eRC")} '
                            f'({scintillator})...', detector.log)
            times = detector.get_lm_data(scintillator, 'SecondsOfDay')
            energies = detector.get_lm_data(scintillator, 'energy')
            wallclock = detector.get_lm_data(scintillator, 'wc')
            count_scints = None

        if detector.processed:
            times = times - detector.first_sec

        # Finding potential events with the search algorithm in find_short_events
        potential_events = find_short_events(detector, modes, scintillator, rollgap, times, energies)

        if len(potential_events) > 0:
            tl.print_logger('Generating scatter plots and event files...', detector.log)

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
            print('', file=detector.log)
            print('Potential short events:', file=detector.log)
            for j in range(len(potential_events)):
                # Stops making plots/event files/log entries after max reached
                # This is to prevent the program getting stuck on lab test days (with hundreds of thousands of "events")
                if (plots_made + plots_already_made) >= params.MAX_PLOTS_PER_SCINT:
                    print(f'Max number of loggable events ({params.MAX_PLOTS_PER_SCINT}) reached.',
                          file=detector.log)
                    break

                if not modes['gui']:
                    print(f'{plots_made}/{max_plots}', end='\r')
                elif modes['gui'] and filecount_switch:
                    print(f'Making {max_plots} plots...')
                    filecount_switch = False

                event = potential_events[j]
                event.number = j + 1 + plots_already_made

                # Logs the event
                start_second = times[event.start] - params.SEC_PER_DAY if (
                        times[event.start] > params.SEC_PER_DAY) else times[event.start]
                print(f'{dt.datetime.utcfromtimestamp(times[event.start] + detector.first_sec)} UTC '
                      f'({start_second} seconds of day) - weather: {tl.weather_from_score(event.weather_subscore)}',
                      file=detector.log)
                print(f'    Rank: {event.rank}, Subscores: [Length: {event.len_subscore}, '
                      f'Clumpiness: {event.clumpiness_subscore}, HEL: {event.hel_subscore}]\n',
                      file=detector.log)

                # Finds the file(s) that the event occurred in
                event.lm_files = find_se_files(detector, event, times, count_scints)

                # Finds the trace(s) associated with the event and aligns them
                event.traces = find_se_traces(detector, event, trace_dict, times, count_scints)

                # Makes the scatter plot for the event
                make_se_scatterplot(detector, event, times, energies, count_scints)

                # Makes the json file for the event
                make_se_json(detector, event, times, energies, wallclock, count_scints)

                plots_made += 1

            if not modes['gui']:
                print(f'{plots_made}/{max_plots}\n', end='\r')

            event_numbers[scintillator] = plots_made + plots_already_made

            tl.print_logger('', detector.log)

        # In combo mode, we only need to run through this loop once
        if modes['combo']:
            break

    return event_numbers


# Combines time data from several scintillators (if applicable) and cuts out counts below a certain energy
def long_event_cutoff(detector, modes, chunk=None):
    # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
    operating_obj = chunk if chunk is not None else detector  # operating_obj should contain data
    # Checks to see that data is present in all preferred scintillators. Otherwise, uses the default scintillator.
    long_event_scintillators = operating_obj.long_event_scint_list
    for scintillator in operating_obj.long_event_scint_list:
        if not operating_obj.data_present_in(scintillator):
            long_event_scintillators = [operating_obj.default_scintillator]
            break

    times = []
    energies = []
    # Checks to see if scintillators have all been calibrated. If one or more aren't, no counts are cut out
    all_calibrated = True
    for scintillator in long_event_scintillators:
        times.append(operating_obj.get_lm_data(scintillator, 'SecondsOfDay'))
        if modes['skcali']:
            energies.append(operating_obj.get_lm_data(scintillator, 'energy'))
        else:
            try:
                energies.append(operating_obj.get_lm_data(scintillator, 'energy', to_mev=True))
            except ValueError:
                all_calibrated = False
                energies.append(operating_obj.get_lm_data(scintillator, 'energy'))

    times = np.concatenate(times)
    energies = np.concatenate(energies)
    # Removes entries that are below a certain cutoff energy
    if not all_calibrated or modes['skcali']:
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
    while low <= high:
        mid = low + (high - low) // 2
        if abs_zscores[sorting_order[mid]] >= 3:
            high = mid - 1
        else:
            low = mid + 1

    # Recalculates mue without outliers
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

        # Fitting curve is linear if either of the windows contains a zero-value bin
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
def find_long_events(modes, day_bins, hist_allday, mue, sigma):
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

    return potential_glows


# Finds the list mode file(s) associated with a long event
def find_le_files(detector, event):
    lm_files = {}
    long_event_scintillators = detector.long_event_scint_list
    for scintillator in detector.long_event_scint_list:
        if not detector.data_present_in(scintillator):
            long_event_scintillators = [detector.default_scintillator]
            break

    for scintillator in long_event_scintillators:
        lm_files[scintillator] = detector.find_lm_file(scintillator, event.start_sec)

    return lm_files


# Makes the histogram subplots for long events
def make_hist_subplot(ax, event, day_bins, hist_allday, mue, sigma):
    left = 0 if (event.peak_index - params.LE_SUBPLOT_PADDING) < 0 else (event.peak_index - params.LE_SUBPLOT_PADDING)
    right = (len(day_bins) - 2) if (event.peak_index + params.LE_SUBPLOT_PADDING) > (len(day_bins) - 2) else \
        (event.peak_index + params.LE_SUBPLOT_PADDING)

    sub_bins = day_bins[left:right]
    ax.bar(sub_bins, hist_allday[left:right], alpha=params.LE_SUBPLOT_BAR_ALPHA,
           color=params.LE_SUBPLOT_BAR_COLOR, width=params.BIN_SIZE)
    ax.set_xlabel('Seconds of Day (UT)')
    ax.set_ylabel('Counts/bin')
    ax.plot(sub_bins, mue[left:right] + params.FLAG_THRESH * sigma[left:right],
            color=params.LE_THRESH_LINE_COLOR, linestyle='dashed', linewidth=2)
    ax.grid(True)


# Runs the long event search
def long_event_search(detector, modes, times, existing_hist=None):
    # Makes one bin for every binsize seconds of the day (plus around 300 seconds more for the next day)
    day_bins = np.arange(0, params.SEC_PER_DAY + 200 + params.BIN_SIZE, params.BIN_SIZE)

    # Creates numerical values for histograms using numpy
    if existing_hist is not None:
        hist_allday = existing_hist
    else:
        hist_allday, _ = np.histogram(times, bins=day_bins)

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

    # Finding potential events with the search algorithm in find_long_events
    potential_glows = find_long_events(modes, day_bins, hist_allday, mue, sigma)

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

    ax1.bar(day_bins[:-1], hist_allday, alpha=params.LE_MAIN_BAR_ALPHA,
            color=params.LE_MAIN_BAR_COLOR, width=params.BIN_SIZE)
    ax1.set_xlabel('Seconds of Day (UT)')
    ax1.set_ylabel('Counts/bin')
    ax1.plot(day_bins[:-1], mue + params.FLAG_THRESH * sigma, color=params.LE_THRESH_LINE_COLOR, linestyle='dashed')

    # Creates legend
    allday_data = mpatches.Patch(color=params.LE_MAIN_BAR_COLOR, label='All Energies')
    allday_thresh_sigma = mpatches.Patch(color=params.LE_THRESH_LINE_COLOR,
                                         label=f'{params.FLAG_THRESH} Sigma Above All Energies', linestyle='dashed')
    ax1.legend(handles=[allday_data, allday_thresh_sigma, ], bbox_to_anchor=(1.05, 1), loc=1,
               borderaxespad=0.)
    ax1.grid(True)

    # Making event files and subplots (if they exist)
    if len(potential_glows) > 0:
        tl.print_logger('Generating event files...', detector.log)
        print('', file=detector.log)
        print('Potential glows:', file=detector.log)

        event_path = f'{detector.get_results_loc()}/event_files/long_events'
        tl.make_path(event_path)

        files_made = 0
        filecount_switch = True
        for i in range(len(potential_glows)):
            if not modes['gui']:
                print(f'{files_made}/{len(potential_glows)}', end='\r')
            elif modes['gui'] and filecount_switch:
                print(f'Making {len(potential_glows)} event files...')
                filecount_switch = False

            glow = potential_glows[i]
            info = (f'{dt.datetime.utcfromtimestamp(glow.start_sec + detector.first_sec)} UTC ({glow.start_sec} '
                    f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                    f'{glow.highest_score}')

            # Logging the event
            print(info, file=detector.log)

            # Making the event file
            event_file = open(f'{event_path}/{detector.date_str}_event{i + 1}_zscore'
                              f'{int(glow.highest_score)}.txt', 'w')
            print(info, file=event_file)
            glow.lm_files = find_le_files(detector, glow)
            for scintillator in glow.lm_files:
                print(f'{scintillator}: {glow.lm_files[scintillator]}', file=event_file)

            event_file.close()
            files_made += 1

        print('', file=detector.log)

        if not modes['gui']:
            print(f'{files_made}/{len(potential_glows)}\n', end='\r')

        # Sorts the glows in descending order depending on their highest z-scores
        potential_glows = sorted(potential_glows, key=lambda x: -x.highest_score)  # Negative for descending order sort

        # Makes the histogram subplots
        ax_list = [ax2, ax3, ax4, ax5]
        for i in range(4):
            try:
                glow = potential_glows[i]
                make_hist_subplot(ax_list[i], glow, day_bins, hist_allday, mue, sigma)
            except IndexError:
                break

    else:
        tl.print_logger(f'There were no potential glows on {detector.full_date_str}', detector.log)

    plt.tight_layout()

    # Saves the histogram(s):
    tl.print_logger('Saving Histogram...', detector.log)
    hist_path = f'{detector.get_results_loc()}'
    tl.make_path(hist_path)
    plt.savefig(f'{hist_path}/{detector.date_str}_histogram.png', dpi=500)
    plt.close(figure)


# A little generator for dividing lists into even-ish parts (used in low memory mode)
def divide_list(lst, num_parts):
    index = int(len(lst) / num_parts)
    part = 1
    if index == 0:
        yield lst
        lst = []
        part += 1

    while part <= num_parts:
        if part < num_parts:
            yield lst[:index]
            lst = lst[index:]
        else:
            yield lst

        part += 1


# Gets necessary info from command line args and then runs the program
def main():
    if len(sys.argv) >= 4:
        first_date = str(sys.argv[1])
        second_date = str(sys.argv[2])
        unit = str(sys.argv[3])
    else:
        print('Please provide a first date, a second date, and a unit name.')
        exit()

    # Makes sure inputs are valid
    if not first_date.isdigit() or not second_date.isdigit() \
            or len(first_date) != 6 or len(second_date) != 6:
        print('Invalid date(s).')
        exit()
    elif int(second_date) < int(first_date):
        print('Not a valid date range.')
        exit()

    if len(sys.argv) > 4:
        mode_info = sys.argv[4:]
    else:
        mode_info = []

    program(first_date, second_date, unit, mode_info)


# Main program function
def program(first_date, second_date, unit, mode_info):
    modes = get_modes(mode_info)

    # Makes a list of all the dates on the requested range
    requested_dates = tl.make_date_list(first_date, second_date)

    # Looping through the dates
    for date_str in requested_dates:
        low_memory_mode = False
        print('')
        print(f'{tl.short_to_full_date(date_str)}:')

        # Initializes the detector object
        print('Importing data...')
        try:
            detector = get_detector(unit, date_str, print_feedback=True)
        except ValueError:
            print('Not a valid detector.')
            exit()

        # Tells the detector to use processed data if the user asks for it
        if modes['processed']:
            detector.use_processed()

        # Tells the detector to use custom import/export directories if the user asks for it
        if modes['custom']:
            index = mode_info.index('-c')
            if index + 2 < len(mode_info):
                import_index = index + 1
                if mode_info[import_index] != 'none':
                    if mode_info[import_index] != '/':
                        detector.set_import_loc(mode_info[import_index])

                export_index = index + 2
                if mode_info[export_index] != 'none':
                    if mode_info[export_index] != '/':
                        detector.set_results_loc(mode_info[export_index])

        # Logs relevant data files and events in a .txt File
        log_path = f'{detector.get_results_loc()}'
        tl.make_path(log_path)
        log = open(f'{log_path}/log.txt', 'w')
        print(f'{tl.short_to_full_date(date_str)}:', file=log)

        # Normal operating mode
        try:
            # Imports the data
            if modes['pickle']:  # In pickle mode: reads the pickle file, or exports one if it doesn't exist yet
                pickle_paths = glob.glob(f'{detector.get_results_loc()}/detector.pickle')
                if len(pickle_paths) > 0:
                    detector = tl.unpickle_detector(pickle_paths[0])
                    detector.log = log
                else:
                    detector.log = log
                    detector.import_data(gui=modes['gui'])
                    tl.pickle_detector(detector, 'detector')
            else:
                detector.log = log
                detector.import_data(gui=modes['gui'])

            # raise MemoryError  # for low memory mode testing

            print('')
            print('Done.')

            # Trace Search
            if not modes['skshort']:
                tl.print_logger('\n', detector.log)
                tl.print_logger('Filtering traces...', detector.log)
                trace_dict = trace_search(detector)
                tl.print_logger('Done.', detector.log)
            else:
                trace_dict = {scintillator: [] for scintillator in detector}

            # Checks to see that necessary list mode data is present
            if not detector.data_present_in(detector.default_scintillator):
                tl.print_logger('\n', detector.log)
                tl.print_logger('No/missing necessary data for specified day.', detector.log)
                raise FileNotFoundError('data missing for one or more scintillators.')

            # Calibrates each scintillator
            if not modes['skcali']:
                tl.print_logger('\n', detector.log)
                tl.print_logger('Calibrating scintillators and generating energy spectra...', detector.log)

                # Calling the calibration algorithm
                detector.calibrate(plot_spectra=True, make_template=modes['template'])

                tl.print_logger('Done.', detector.log)

            # Short event search
            if not modes['skshort']:
                tl.print_logger('\n', detector.log)
                tl.print_logger('Starting search for short events...', detector.log)

                # Calling the short event search algorithm
                short_event_search(detector, modes, trace_dict)

                tl.print_logger('Done.', detector.log)

            # Long event search
            if not modes['skglow']:
                tl.print_logger('\n', detector.log)
                tl.print_logger('Starting search for glows...', detector.log)
                # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
                le_times = long_event_cutoff(detector, modes)

                # Calling the long event search algorithm
                long_event_search(detector, modes, le_times)

                tl.print_logger('Done.', detector.log)

        except MemoryError:
            low_memory_mode = True

        except FileNotFoundError:  # Missing necessary data
            pass

        except Exception as ex:  # Logging errors
            tl.print_logger('\n', log)
            tl.print_logger(f'Search could not be completed due to the following error: {ex}', log)
            tl.print_logger('See error log for details.', log)
            with open(f'{log_path}/err.txt', 'w') as err_file:
                err_file.write(traceback.format_exc())

        # Low memory mode
        if low_memory_mode:
            try:
                tl.print_logger('', detector.log)
                tl.print_logger('Not enough memory. Entering low memory mode...', detector.log)
                # Measures the total combined size of all the data files
                total_file_size = detector.get_fileset_size()
                # Clears leftover data (just to be sure)
                detector.clear(clear_filelists=False)

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
                    chunk = get_detector(unit, date_str, print_feedback=True)
                    chunk.log = log

                    if modes['processed']:
                        chunk.use_processed()

                    if modes['custom']:
                        index = mode_info.index('-c')
                        if index + 2 < len(mode_info):
                            import_index = index + 1
                            if mode_info[import_index] != 'none':
                                if mode_info[import_index] != '/':
                                    chunk.set_import_loc(mode_info[import_index])

                            export_index = index + 2
                            if mode_info[export_index] != 'none':
                                if mode_info[export_index] != '/':
                                    chunk.set_results_loc(mode_info[export_index])

                    chunk_list.append(chunk)

                chunk_scint_list = chunk_list[0].scint_list

                for scintillator in detector:
                    lm_filelist_chunks = divide_list(detector.get_attribute(scintillator, 'lm_filelist'),
                                                     num_chunks)
                    trace_filelist_chunks = divide_list(detector.get_attribute(scintillator, 'trace_filelist'),
                                                        num_chunks)
                    for chunk in chunk_list:
                        chunk.set_attribute(scintillator, 'lm_filelist', next(lm_filelist_chunks), deepcopy=False)
                        chunk.set_attribute(scintillator, 'trace_filelist', next(trace_filelist_chunks), deepcopy=False)

                # Imports data to each chunk and then pickles the chunks (and checks that data is actually present)
                print('Importing data...')

                chunk_num = 1
                chunk_path_list = []
                has_data = True
                # Temporary pickle feature for low memory mode. REMOVE WHEN PROGRAM IS FINISHED
                pickled_chunk_paths = glob.glob(f'{detector.get_results_loc()}/Results/{unit}/{date_str}/chunk*.pickle')
                pickled_chunk_paths.sort()
                if modes['pickle'] and len(pickled_chunk_paths) > 0:
                    chunk_path_list = pickled_chunk_paths
                else:
                    # Keeps timings consistent between chunks
                    passtime_dict = {scintillator: chunk_list[0].get_attribute(scintillator, 'passtime', deepcopy=True)
                                     for scintillator in chunk_scint_list}

                    for chunk in chunk_list:
                        # Updates chunk to include previous chunk's passtime
                        if chunk != chunk_list[0]:
                            for scintillator in chunk:
                                chunk.set_attribute(scintillator, 'passtime', passtime_dict[scintillator],
                                                    deepcopy=True)

                        tl.print_logger('', detector.log)
                        tl.print_logger(f'Chunk {chunk_num} (of {num_chunks}):', detector.log)
                        chunk.import_data(existing_filelists=True, gui=modes['gui'])

                        # Checking that data is present in the necessary scintillators
                        if not chunk.data_present_in(chunk.default_scintillator):
                            has_data = False

                        # Makes a full list of filetime extrema for long event search
                        # Also updates passtime_dict for the next chunk
                        for scintillator in chunk:
                            extrema = detector.get_attribute(scintillator, 'lm_file_ranges')
                            extrema += chunk.get_attribute(scintillator, 'lm_file_ranges')
                            detector.set_attribute(scintillator, 'lm_file_ranges', extrema)
                            if chunk != chunk_list[-1]:
                                passtime_dict[scintillator] = chunk.get_attribute(scintillator, 'passtime',
                                                                                  deepcopy=True)

                        # Pickle chunk and add its path to the list
                        chunk_path_list.append(tl.pickle_chunk(chunk, f'chunk{chunk_num}'))

                        # Eliminates the chunk from active memory
                        del chunk
                        gc.collect()
                        chunk_list[chunk_num - 1] = 0

                        chunk_num += 1

                print('')
                print('Done.')

                # Trace Search
                chunk_trace_dicts = {}
                if not modes['skshort']:
                    tl.print_logger('\n', detector.log)
                    tl.print_logger('Filtering traces...', detector.log)
                    for chunk_path in chunk_path_list:
                        chunk = tl.unpickle_chunk(chunk_path)
                        trace_dict = trace_search(chunk)
                        chunk_trace_dicts[chunk_path] = trace_dict

                    tl.print_logger('Done.', detector.log)
                else:
                    for chunk_path in chunk_path_list:
                        chunk_trace_dicts[chunk_path] = {scintillator: [] for scintillator in chunk_scint_list}

                # Skips the day if necessary list mode data is missing in any of the chunks
                if not has_data:
                    tl.print_logger('No/missing necessary data for specified day.', detector.log)
                    raise FileNotFoundError('data missing for one or more scintillators.')

                # Calibrates each scintillator
                if not modes['skcali']:
                    tl.print_logger('\n', detector.log)
                    tl.print_logger('Calibrating scintillators and generating energy spectra...', detector.log)
                    # Keeps track of spectra contributions for each chunk
                    existing_spectra = {scintillator: np.array([]) for scintillator in chunk_scint_list}
                    energy_bins = np.arange(0.0, detector.calibration_params['bin_range'],
                                            detector.calibration_params['bin_size'])
                    for chunk_path in chunk_path_list:
                        chunk = tl.unpickle_chunk(chunk_path)
                        for scintillator in chunk_scint_list:
                            energies = chunk.get_lm_data(scintillator, 'energy')
                            chunk_hist, _ = np.histogram(energies, bins=energy_bins)
                            if len(existing_spectra[scintillator]) == 0:
                                existing_spectra[scintillator] = chunk_hist
                            else:
                                existing_spectra[scintillator] = existing_spectra[scintillator] + chunk_hist

                        del chunk

                    # Calling the calibration algorithm
                    detector.calibrate(existing_spectra=existing_spectra, plot_spectra=True,
                                       make_template=modes['template'])

                    tl.print_logger('Done.', detector.log)

                # Short event search
                if not modes['skshort']:
                    tl.print_logger('\n', detector.log)
                    tl.print_logger('Starting search for short events...', detector.log)
                    tl.print_logger('Warning! In low memory mode, short events '
                                    'will be ranked on a per-chunk basis.',
                                    detector.log)
                    tl.print_logger('', detector.log)
                    prev_event_numbers = {scintillator: 0 for scintillator in chunk_scint_list}

                    chunk_num = 1
                    for chunk_path in chunk_path_list:
                        chunk = tl.unpickle_chunk(chunk_path)
                        chunk.log = log
                        tl.print_logger(f'Chunk {chunk_num} (of {num_chunks}):', detector.log)

                        # Calling the short event search algorithm
                        prev_event_numbers = short_event_search(chunk, modes,
                                                                chunk_trace_dicts[chunk_path],
                                                                prev_event_numbers=prev_event_numbers)

                        del chunk
                        gc.collect()
                        chunk_num += 1

                    tl.print_logger('Done.', detector.log)

                # Long event search
                if not modes['skglow']:
                    tl.print_logger('\n', detector.log)
                    tl.print_logger('Starting search for glows...', detector.log)
                    le_hist = np.array([])
                    day_bins = np.arange(0, params.SEC_PER_DAY + 200 + params.BIN_SIZE, params.BIN_SIZE)
                    for chunk_path in chunk_path_list:
                        chunk = tl.unpickle_chunk(chunk_path)
                        # Converts energy channels to MeV using the locations of peaks/edges
                        # obtained during calibration
                        le_times = long_event_cutoff(detector, modes, chunk)

                        # Histograms the counts from each chunk and combines them with the main one
                        chunk_hist, _ = np.histogram(le_times, bins=day_bins)
                        le_hist = chunk_hist if len(le_hist) == 0 else le_hist + chunk_hist

                        del chunk
                        gc.collect()

                    # Calling the long event search algorithm
                    long_event_search(detector, modes, np.array([]), existing_hist=le_hist)

                    tl.print_logger('Done.', detector.log)

                log.close()
                # Deletes chunk .pickle files
                # REMOVE CONDITIONAL STATEMENT WHEN PROGRAM IS DONE
                if not modes['pickle']:
                    for chunk_path in chunk_path_list:
                        os.remove(chunk_path)

            except MemoryError:
                tl.print_logger('\n', log)
                tl.print_logger('Cannot complete search. Too little memory available on system.', log)

            except FileNotFoundError:  # Missing necessary data
                pass

            except Exception as ex:  # Logging errors
                tl.print_logger('\n', log)
                tl.print_logger(f'Search could not be completed due to the following error: {ex}', log)
                tl.print_logger('See error log for details.', log)
                with open(f'{log_path}/err.txt', 'w') as err_file:
                    err_file.write(traceback.format_exc())

        del detector
        log.close()
        gc.collect()


if __name__ == '__main__':
    main()
