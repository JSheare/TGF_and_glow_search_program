import sys
import warnings
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import search_classes as sc
import search_module as sm
import pandas as pd

# Disables all warnings in the terminal because they're mostly just annoying
warnings.filterwarnings('ignore')

# Establish dates from inputs 1 & 2
input1 = str(sys.argv[1])
year_int = int('20' + input1[:2])
month_int = int(input1[2:4])
day_int = int(input1[-2:])
year_month = input1[:4]

input2 = str(sys.argv[2])
year_int2 = int('20' + input2[:2])
month_int2 = int(input2[2:4])
day_int2 = int(input2[-2:])
year_month2 = input2[:4]

# Input 3 is where the user specifies which set of detector data is desired
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

# Makes a dictionary of all the requested years, the requested months in each year,
# and the requested days in each month
requested_dates = {}
current_year = year_int
# For the years
while True:
    if current_year == year_int2:
        requested_dates.update({str(current_year): []})
        break
    else:
        requested_dates.update({str(current_year): []})
        current_year += 1

current_month = month_int
# For the months and days
for i in requested_dates:
    requested_months = []
    # For requested dates that are all in the same month
    if year_int == year_int2 and month_int == month_int2:
        requested_months.append([current_month, day_int, day_int2])
    # For everything else
    else:
        while True:
            if i == str(year_int) and current_month == month_int:
                requested_months.append([current_month, day_int, sm.days_per_month(current_month, int(i))])
                current_month += 1
            elif i == str(year_int2) and current_month == month_int2:
                requested_months.append([current_month, 1, day_int2])
                break
            elif current_month > 12:
                current_month = 1
                break
            else:
                requested_months.append([current_month, 1, sm.days_per_month(current_month, int(i))])
                current_month += 1

    current_month = 1
    requested_dates.update({i: requested_months})

for year in requested_dates:  # Loops over all requested years
    requested_months = requested_dates[year]
    YEAR = int(year)
    YEAR_str = year
    for month in requested_months:  # Loops over each requested month of the given year
        MONTH = month[0]
        MONTH_str = f'0{month[0]}' if len(str(month[0])) < 2 else str(month[0])
        day1 = month[1]
        day2 = month[2]
        for day in range(day1, day2+1):  # Loops over each requested day of the given month
            DAY = day
            DAY_str = f'0{day}' if len(str(day)) < 2 else str(day)
            full_day_string = f'{YEAR_str[2:4]}{MONTH_str}{DAY_str}'  # In format yymmdd
            date_timestamp = dt.date(YEAR, MONTH, DAY)  # In format yyyy-mm-dd
            print(f'\n{date_timestamp}:')

            # Logs relevant data files and events in a .txt File
            log_path = f'{sm.results_loc()}Results/{unit}/{full_day_string}/'
            sm.path_maker(log_path)
            log = open(f'{log_path}log.txt', 'w')

            # EPOCH time conversions
            first_sec = (dt.datetime(YEAR, MONTH, DAY, 0, 0) - dt.datetime(1970, 1, 1)).total_seconds()
            sm.print_logger(f'The first second of {YEAR}-{MONTH}-{DAY} is: {int(first_sec)}', log)
            sm.print_logger(f'The last second of {YEAR}-{MONTH}-{DAY} is: {int(first_sec + 86400)}', log)

            # Initializes the detector object and imports the data
            print('Importing data...')
            detector = sc.Detector(unit, first_sec, log, modes)
            detector.data_importer()

            if detector.GODOT or detector.THOR:
                if len(detector.scintillators['NaI']['filelist']) == 0 or\
                        len(detector.scintillators['LP']['filelist']) == 0:
                    print('\n\n')
                    print('\n', file=detector.log)
                    sm.print_logger('No/Missing data for specified day.', detector.log)
                    print('\n')
                    continue
            else:  # Santis instrument
                if len(detector.scintillators['LP']['filelist']) == 0:
                    print('\n\n')
                    print('\n', file=detector.log)
                    sm.print_logger('No/Missing data for specified day.', detector.log)
                    print('\n')
                    continue

            print('\n')
            print('Done.')

            # Calibrates each scintillator
            lp_channels = []
            nai_channels = []
            if not skcali:
                print('\n')
                sm.print_logger('Calibrating Scintillators and generating energy spectra...', detector.log)
                lp_channels, nai_channels = detector.spectra_maker()
                sm.print_logger('Done.', detector.log)

    # Short event algorithm starts here:
            if not skshort:
                sm.print_logger('\n', detector.log)
                sm.print_logger('Starting search for short events...', detector.log)
                sm.print_logger('\n', detector.log)

                # Parameters:
                rollgap = 4 if aircraft else 4
                event_time_spacing = 1e-3  # 1 millisecond
                event_min_counts = 10
                noise_cutoff_energy = 300
                min_noise_counts = 3
                channel_range_width = 300
                channel_ratio = 0.5

                for scintillator in detector.scintillators:
                    if scintillator != 'LP' and not detector.plastics:
                        continue
                    elif scintillator == 'NaI' and detector.plastics:
                        continue

                    sm.print_logger(f'For eRC {detector.scintillators[scintillator]["eRC"]} '
                                    f'({scintillator}):', detector.log)
                    times = detector.attribute_retriever(scintillator, 'time')
                    energies = detector.attribute_retriever(scintillator, 'energy')
                    filelist = detector.attribute_retriever(scintillator, 'filelist')
                    filetime_extrema = detector.attribute_retriever(scintillator, 'filetime_extrema')

                    # Checks for an event by looking for a certain number of counts (rollgap + 1) in a small timeframe
                    potential_event_list = []
                    event_start = 0
                    event_length = 0
                    event_time = 0
                    total_potential_events = 0
                    total_threshold_reached = 0

                    rolled_array = np.roll(times, -rollgap)
                    interval = rolled_array - times
                    event_start_time = 0
                    for i in range(len(interval)):
                        if 0 < interval[i] <= event_time_spacing:
                            # Records the beginning index of a potential event
                            if event_length == 0:
                                event_start = i
                                event_start_time = times[i]
                                event_length = 1 + rollgap  # 1 for first count, rollgap for the others
                                event_time = rolled_array[i] - event_start_time
                                total_potential_events += 1
                                print(f'Potential event (#{total_potential_events})')
                            # Measures the length of a potential event
                            else:
                                event_length += 1
                                event_time = rolled_array[i] - event_start_time

                            # Counts the total number of times that the detection threshold was reached
                            total_threshold_reached += 1

                        # Records the rough length of a potential event
                        if interval[i] > event_time_spacing and event_length > 0:
                            # Keeps potential event if it's longer than the specified minimum number of counts
                            if event_length >= event_min_counts:
                                print(f'Potential event length: {event_time} seconds, {event_length} counts')
                                potential_event = sc.ShortEvent(event_start, event_length, scintillator)
                                potential_event_list.append(potential_event)

                            if event_length < event_min_counts:
                                print(f'Potential event removed due to insufficient length ({event_length} counts)')

                            event_start = 0
                            event_length = 0
                            event_time = 0

                    # Eliminates noisy events
                    f_potential_event_list = []
                    for event in potential_event_list:
                        event_energies = energies[event.start:event.stop]
                        event_length = len(event_energies)

                        # Checks that there are fewer counts in the higher energy channels than the low energy channels
                        # High/low channel ratio is not checked for THOR
                        good_channel_ratio = True if detector.THOR else False
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
                        # High/low channel ratio is not checked for events with < 30 counts
                        else:
                            good_channel_ratio = True

                        # Eliminates the low energy events, which are likely just noise
                        # These 5 lines can probably be condensed into one, but I'm not going to do it right now
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
                                print('Potential event removed due to noise (bad high-to-low channel ratio)')
                            elif not is_greater_than_thresh and good_channel_ratio:
                                print('Potential event removed due to noise (minimum energy threshold not reached)')
                            else:
                                print('Potential event removed due to noise')

                    potential_event_list = f_potential_event_list

                    sm.print_logger(f'\n{len(potential_event_list)} potential events recorded', detector.log)
                    sm.print_logger(f'Detection threshold reached {total_threshold_reached} times', detector.log)
                    print('\n', file=detector.log)

                    if len(potential_event_list) > 0:
                        print('\n')
                        sm.print_logger('Generating scatter plots and event files...', detector.log)
                        print('\n', file=detector.log)
                        print('Potential short events:', file=detector.log)
                        for event in potential_event_list:
                            start_second = times[event.start] - 86400 if times[event.start] > 86400 else \
                                times[event.start]
                            print(f'{dt.datetime.utcfromtimestamp(times[event.start] + first_sec)} UTC '
                                  f'({start_second} seconds of day)', file=detector.log)

                        print('\n', file=detector.log)

                        # Makes scatter plots of the resulting potential events
                        # Subplot timescales
                        ts1 = 1e-4  # 100 microseconds
                        ts2 = 0.005  # 5 milliseconds
                        ts3 = 0.1  # 100 milliseconds
                        ts_list = [ts1, ts2, ts3]

                        plots_made = 0
                        for i in range(len(potential_event_list)):
                            print(f'{plots_made}/{len(f_potential_event_list)}', end='\r')
                            event = potential_event_list[i]
                            filelist, filetime_extrema = event.scatterplot_maker(ts_list, detector, i+1,
                                                                                 filelist, filetime_extrema)
                            plots_made += 1
                            print(f'{plots_made}/{len(potential_event_list)}', end='\r')

                        sm.print_logger('Done.', detector.log)
                        sm.print_logger('\n', detector.log)
                    else:
                        print('\n')

    # Glow search algorithm starts here
            if not skglow:
                sm.print_logger('Starting search for glows...', detector.log)

                # Converts energy channels to MeV using the locations of peaks/edges obtained during calibration
                # This code could probably be cleaned up
                if detector.good_lp_calibration and not skcali:
                    lp_times = detector.attribute_retriever('LP', 'time')
                    lp_energies = sm.channel_to_mev(detector.attribute_retriever('LP', 'energy'), lp_channels, 'LP')
                else:
                    lp_times = detector.attribute_retriever('LP', 'time')
                    lp_energies = detector.attribute_retriever('LP', 'energy')

                if detector.SANTIS:
                    nai_times = np.array([])
                    nai_energies = np.array([])
                elif not skcali:
                    nai_times = detector.attribute_retriever('NaI', 'time')
                    nai_energies = sm.channel_to_mev(detector.attribute_retriever('NaI', 'energy'), nai_channels, 'NaI')
                else:
                    nai_times = detector.attribute_retriever('NaI', 'time')
                    nai_energies = detector.attribute_retriever('NaI', 'energy')

                # Combines large plastic and NaI data for GODOT
                scint_list = []  # This list should really be made into a detector attribute
                if detector.GODOT:
                    scint_list = ['NaI', 'LP']
                    times = np.append(lp_times, nai_times)
                    energies = np.append(lp_energies, nai_energies)
                elif detector.THOR:
                    scint_list = ['NaI']
                    times = nai_times
                    energies = nai_energies
                else:
                    scint_list = ['LP']
                    times = detector.attribute_retriever('LP', 'time')
                    energies = detector.attribute_retriever('LP', 'energy')

                # Removes entries that are below a certain cutoff energy
                if (not detector.good_lp_calibration and (detector.GODOT or detector.SANTIS)) or skcali:
                    print('Potentially inaccurate large plastic calibration, beware radon washout!', file=detector.log)
                else:
                    energy_cutoff = 1.9  # MeV
                    cut_indices = np.where(energies < energy_cutoff)
                    times = np.delete(times, cut_indices)
                    energies = np.delete(energies, cut_indices)

                if detector.processed:
                    times = times - first_sec

                # Makes one bin for every binsize seconds of the day (plus around 300 seconds more for the next day)
                binsize = 4 if aircraft else 10
                day_bins = np.arange(0, 86700 + binsize, binsize)

                # Creates numerical values for histograms using numpy
                hist_allday, bins_allday = np.histogram(times, bins=day_bins)

                # Calculates mean and z-scores
                hist_allday_nz = hist_allday[hist_allday.nonzero()]
                mue = sum(hist_allday_nz) / len(hist_allday_nz)  # Mean number of counts/bin
                sigma = np.sqrt(mue)  # Standard deviation of the distribution

                # Removes outliers (3 sigma) and recalculates mue
                abs_zscores = np.abs((hist_allday_nz - mue) / sigma)
                hist_allday_nz = hist_allday_nz[np.argsort(abs_zscores)]
                hist_allday_nz = hist_allday_nz[::-1]
                zscores = np.sort(abs_zscores)[::-1]
                for i in range(len(abs_zscores)):
                    if abs_zscores[i] < 3:
                        hist_allday_nz = hist_allday_nz[i:]
                        break

                mue = sum(hist_allday_nz) / len(hist_allday_nz)
                sigma = np.sqrt(mue)

                z_scores = np.array([])  # The z-scores themselves
                z_flags = np.array([])
                p = 0
                # Flags only those z-scores > 5
                for i in range(len(hist_allday)):
                    if hist_allday[i] > mue:  # Peak
                        z_scores = np.append(z_scores, ((hist_allday[i] - mue) / sigma))
                        if z_scores[i] >= 5:
                            z_flags = np.append(z_flags, i)
                            p += 1
                    elif hist_allday[i] < mue:  # Valley
                        z_scores = np.append(z_scores, ((hist_allday[i] - mue) / sigma))

                # Redefines z_flag to only contain flags from the start to p-1
                # yeah, but why?
                z_flags = z_flags[0:p].astype(int)

                # Sorts z-flags into actual potential glows
                glow_start = 0
                glow_length = 0
                glow_time = 0
                previous_time = 0
                potential_glow_list = []
                for flag in z_flags:
                    if glow_length == 0:
                        glow_start = flag
                        glow_length += 1
                    if glow_length > 0 and day_bins[flag] - binsize == previous_time:
                        glow_length += 1
                    if glow_length > 0 and day_bins[flag] - binsize > previous_time:
                        # Makes glow object and fills it out
                        glow = sc.PotentialGlow(glow_start, glow_length)
                        glow.highest_score = glow.highest_zscore(z_scores)
                        glow.start_sec, glow.stop_sec = glow.beginning_and_end_seconds(day_bins, binsize)
                        potential_glow_list.append(glow)
                        glow_start = flag
                        glow_length = 1
                    if flag == z_flags[-1]:
                        glow = sc.PotentialGlow(glow_start, glow_length)
                        glow.highest_score = glow.highest_zscore(z_scores)
                        glow.start_sec, glow.stop_sec = glow.beginning_and_end_seconds(day_bins, binsize)
                        potential_glow_list.append(glow)

                    previous_time = day_bins[flag]

                # Checks the bins surrounding a long event while in aircraft mode
                if aircraft:
                    # Number of neighbors on each side of the long event edges to check
                    check_neighbor = 10
                    # Acceptable difference between peak zscore of event and average zscore of surrounding events
                    tolerance = 1

                    num_bins = len(z_scores)
                    f_potential_glow_list = []
                    for glow in potential_glow_list:
                        z_score_sum = 0
                        bins_summed = 0
                        for i in range(check_neighbor):
                            if not glow.start - (i+1) < 0:
                                z_score_sum += z_scores[glow.start - (i+1)]
                                bins_summed += 1

                            if not glow.start + (i+1) > (num_bins - 1):
                                z_score_sum += z_scores[glow.stop + (i+1)]
                                bins_summed += 1

                        score_average = z_score_sum/bins_summed
                        if (np.abs(z_scores[glow.peak_index] - score_average)) >= tolerance:
                            f_potential_glow_list.append(glow)

                    potential_glow_list = f_potential_glow_list

                sm.print_logger('Done.', detector.log)
                if len(potential_glow_list) == 0:
                    sm.print_logger(f'There were no potential glows for the date {date_timestamp}', detector.log)
                else:
                    # Logs potential glows and sorts them in descending order depending on their highest z-score
                    sm.print_logger('\n', detector.log)
                    sm.print_logger('Generating event files...', detector.log)
                    highest_scores = []
                    print('\n', file=detector.log)
                    print('Potential glows:', file=detector.log)

                    eventpath = f'{sm.results_loc()}Results/{detector.unit}/{detector.full_day_string}/event files/' \
                                f'long events/'
                    sm.path_maker(eventpath)
                    event_number = 1
                    files_made = 0
                    for i in range(len(potential_glow_list)):
                        print(f'{files_made}/{len(potential_glow_list)}', end='\r')
                        glow = potential_glow_list[i]
                        highest_score = glow.highest_score
                        highest_scores.append(highest_score)
                        event_frame = pd.DataFrame()
                        print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + first_sec)} UTC ({glow.start_sec} '
                              f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                              f'{highest_score}', file=detector.log)

                        event_file = open(f'{eventpath}{detector.full_day_string}_event{event_number}_zscore'
                                          f'{int(highest_score)}.txt', 'w')
                        print(f'{dt.datetime.utcfromtimestamp(glow.start_sec + first_sec)} UTC ({glow.start_sec} '
                              f'seconds of day), {glow.stop_sec - glow.start_sec} seconds long, highest z-score: '
                              f'{highest_score}', file=event_file)
                        for scint in scint_list:
                            print(f'{scint}:', file=event_file)
                            filelist = detector.attribute_retriever(scint, 'filelist')
                            filetime_extrema = detector.attribute_retriever(scint, 'filetime_extrema')
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

                    print('\n', file=detector.log)
                    sm.print_logger('Done.', detector.log)

                    glow_sorting_order = np.argsort(highest_scores)
                    glow_sorting_order = glow_sorting_order[::-1]
                    potential_glow_list = [potential_glow_list[s] for s in glow_sorting_order]

                    # Establishes detector location based on year and detector name
                    location = sm.location(unit, YEAR)

                    # Plotting the histograms
                    sm.print_logger('\n', detector.log)
                    sm.print_logger('Generating Histogram...', detector.log)
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
                    n, bins, patches = ax1.hist(times, bins=day_bins, alpha=0.5, color='r')
                    ax1.set_xlabel('Seconds of Day (UT)')
                    ax1.set_ylabel('Counts/10-second bin')
                    ax1.axhline(y=(mue + 5 * sigma), color='blue', linestyle='dashed', linewidth=2)

                    # Creates legend
                    allday_data = mpatches.Patch(color='r', label='All Energies')
                    allday_5sigma = mpatches.Patch(color='blue', label='5 Sigma Above All Energies', linestyle='dashed')
                    ax1.legend(handles=[allday_data, allday_5sigma, ], bbox_to_anchor=(1.05, 1), loc=1,
                               borderaxespad=0.)
                    ax1.grid(True)

                    for i in range(4):
                        try:
                            glow = potential_glow_list[i]
                            sm.hist_subplotter(ax_list[i], glow, times, day_bins, mue, sigma)
                        except IndexError:
                            continue

                    plt.tight_layout()
                    plt.plot()

                    # Saves the histograms:
                    hist_path = f'{sm.results_loc()}Results/{unit}/{full_day_string}/'
                    sm.path_maker(hist_path)
                    plt.savefig(f'{hist_path}{full_day_string}_histogram.png', dpi=500)
                    plt.close(figu)
                    sm.print_logger('Done', detector.log)
                    sm.print_logger('\n', detector.log)

            log.close()
