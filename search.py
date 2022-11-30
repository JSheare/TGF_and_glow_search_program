import sys
import warnings
import datetime
import calendar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import search_classes as sc
import search_module as sm

# IMPORTANT: when you want to run this program, use the following command:
# 'python3 search.py yymmdd1 yymmdd2 detector mode'
# ex: 'python3 search.py 210927 210930 GODOT processed' for 2021 Sept 27 to Sept 30 processed GODOT data

# Note: All import/export paths can be changed in search_module.py

# Disables all warnings in the terminal because they're mostly just annoying
warnings.filterwarnings('ignore')

# Establish dates from inputs 1 & 2
input1 = str(sys.argv[1])
year_int = int('20' + input1[:2])
month_int = int(input1[2:4])
day_int = int(input1[-2:])
year_month = input1[:4]

input2 = str(sys.argv[2])
year_int2 = int('20' + input2[:2])  # This is probably redundant, but I'll leave it for now
month_int2 = int(input2[2:4])
day_int2 = int(input2[-2:])
year_month2 = input2[:4]

# Input 3 is where the user specifies which set of detector data is desired
unit = str(sys.argv[3])

try:
    mode = sys.argv[4:]
except IndexError:
    mode = []
    pass

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
    for i in requested_months:  # Loops over each requested month of the given year
        MONTH = i[0]
        MONTH_str = f'0{i[0]}' if len(str(i[0])) < 2 else str(i[0])
        day1 = i[1]
        day2 = i[2]
        for day in range(day1, day2+1):  # Loops over each requested day of the given month
            DAY = day
            DAY_str = f'0{day}' if len(str(day)) < 2 else str(day)
            full_day_string = f'{YEAR_str[2:4]}{MONTH_str}{DAY_str}'  # In format yymmdd
            date_timestamp = datetime.date(YEAR, MONTH, DAY)  # In format yyyy-mm-dd
            print(f'\n{date_timestamp}:')

            # Logs relevant data files and events in a .txt File
            log_path = f'{sm.results_loc()}Results/{unit}/{full_day_string}/'
            sm.path_maker(log_path)
            datetime_logs = open(f'{log_path}log.txt', 'w')

            # Finds the EPOCH time of the first & last second of the day
            first_sec = calendar.timegm((YEAR, MONTH, DAY, 0, 0, 0, -1, -1, -1))
            last_sec = first_sec + 86400
            sm.print_logger(f'The first second of {YEAR}-{MONTH}-{DAY} is: {first_sec}', datetime_logs)
            sm.print_logger(f'The last second of {YEAR}-{MONTH}-{DAY} is: {last_sec}', datetime_logs)

            # Imports the data
            print('Importing data...')
            detector = sc.Detector(unit, full_day_string, mode)
            detector.data_importer(datetime_logs)
            if len(detector.scintillators['NaI']['filelist']) == 0 or (detector.scintillators['LP']['filelist']) == 0:
                print('\n\n')
                print('\n', file=datetime_logs)
                sm.print_logger('No/Missing data for specified day.', datetime_logs)
                print('\n')
                continue

            print('\n')
            print('Done.')

            # Calibrates each scintillator
            print('\n')
            sm.print_logger('Calibrating Scintillators and generating energy spectra graphs...', datetime_logs)
            lp_channels, nai_channels = detector.spectra_maker(date_timestamp, full_day_string, datetime_logs)
            sm.print_logger('Done.', datetime_logs)

    # Short event algorithm starts here:
            sm.print_logger('\n', datetime_logs)
            sm.print_logger('Starting search for short events...', datetime_logs)
            sm.print_logger('\n', datetime_logs)

            # Parameters:
            rollgap = 4
            event_time_spacing = 1e-3  # 1 millisecond
            event_min_counts = 6
            noise_cutoff_energy = 300
            min_noise_counts = 3
            channel_range_width = 300
            channel_ratio = 0.5

            for g in detector.scintillators:
                if g != 'LP' and not detector.plastics:
                    continue
                elif g == 'NaI' and detector.plastics:
                    continue

                sm.print_logger(f'For eRC {detector.scintillators[g]["eRC"]} ({g}):', datetime_logs)
                times = detector.attribute_retriever(g, 'time')
                energies = detector.attribute_retriever(g, 'energy')
                filelist = detector.attribute_retriever(g, 'filelist')

                # Checks for an event by looking for a certain number of counts (rollgap + 1) in a small timeframe
                potential_event_list = []
                event_start = 0
                event_length = 0
                event_time = 0
                total_potential_events = 0
                total_threshold_reached = 0

                interval = np.abs(times - np.roll(times, rollgap))

                # -1's are meant to band-aid a bug. I pretty much know what's causing it but too lazy to fix right now
                for v in range(len(interval)):
                    # Records the beginning index of a potential event
                    if interval[v] < event_time_spacing and event_length == 0:
                        event_start = v
                        event_length += 1
                        event_time += interval[v]
                        total_potential_events += 1
                        print(f'Potential event (#{total_potential_events})')
                    # Measures the length of a potential event
                    if interval[v] < event_time_spacing and event_length >= 1:
                        event_length += 1
                        event_time += interval[v]

                    # Records the rough length of a potential event
                    if interval[v] > event_time_spacing and event_length > 0:
                        # Keeps potential event if it is longer than the specified minimum number of counts
                        if (event_length - 1) >= event_min_counts:
                            print(f'Potential event length: {event_time} seconds, {event_length - 1} counts')
                            potential_event = sc.ShortEvent(event_start, event_length - 1, g)
                            potential_event_list.append(potential_event)

                        if (event_length - 1) < event_min_counts:
                            print(f'Potential event removed due to insufficient length ({event_length - 1} counts)')

                        event_start = 0
                        event_length = 0
                        event_time = 0
                    # Counts the total number of times that the detection threshold was reached
                    if interval[v] < event_time_spacing:
                        total_threshold_reached += 1

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
                        for u in event_energies:
                            if 200 <= u <= (200 + channel_range_width):
                                low_channel_counts += 1

                            if (300 + channel_range_width) <= u <= (300 + 2 * channel_range_width):
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
                    for x in max_indices:
                        max_energies = np.append(max_energies, event_energies[x])

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

                sm.print_logger(f'\n{len(f_potential_event_list)} potential events recorded', datetime_logs)
                sm.print_logger(f'Detection threshold reached {total_threshold_reached} times', datetime_logs)
                print('\n', file=datetime_logs)

                if len(f_potential_event_list) > 0:
                    print('Potential short events:', file=datetime_logs)
                    for event in f_potential_event_list:
                        print(f'{datetime.datetime.utcfromtimestamp(times[event.start] + first_sec)} UTC '
                              f'({times[event.start]} seconds of day)', file=datetime_logs)

                    print('\n', file=datetime_logs)

                    # Makes scatter plots of the resulting potential events
                    print('\n')
                    sm.print_logger('Generating scatter plots...', datetime_logs)

                    # Subplot timescales
                    ts1 = 1e-4  # 100 microseconds
                    ts2 = 0.005  # 5 milliseconds
                    ts3 = 0.1  # 100 milliseconds
                    ts_list = [ts1, ts2, ts3]

                    plots_made = 0
                    for y in range(len(f_potential_event_list)):
                        print(f'{plots_made}/{len(f_potential_event_list)}', end='\r')
                        event = f_potential_event_list[y]
                        new_filelist = event.scatterplot_maker(ts_list, filelist, times, energies,
                                                               y+1, date_timestamp, unit, mode)
                        plots_made += 1
                        print(f'{plots_made}/{len(f_potential_event_list)}', end='\r')

                    sm.print_logger('Done.', datetime_logs)
                    sm.print_logger('\n', datetime_logs)
                    print('\n')
                else:
                    print('\n')

    # Glow search algorithm starts here
            sm.print_logger('Starting search for glows...', datetime_logs)

            # Converts energy channels from volts to MeV using the locations of peaks/edges obtained during calibration
            lp_times = detector.attribute_retriever('LP', 'time')
            lp_energies = sm.channel_to_mev(detector.attribute_retriever('LP', 'energy'), lp_channels, 'LP')
            if detector.SANTIS:
                nai_times = []
                nai_energies = []
            else:
                nai_times = detector.attribute_retriever('NaI', 'time')
                nai_energies = sm.channel_to_mev(detector.attribute_retriever('NaI', 'energy'), nai_channels, 'NaI')

            # Combines large plastic and NaI data for GODOT
            if detector.GODOT:
                times = np.append(lp_times, nai_times)
                energies = np.append(lp_energies, nai_energies)
            elif detector.THOR:
                times = nai_times
                energies = nai_energies
            else:
                times = detector.attribute_retriever('LP', 'time')
                energies = detector.attribute_retriever('LP', 'energy')

            # Removes entries that are below a certain cutoff energy
            if detector.good_lp_calibration:
                energy_cutoff = 1.9  # MeV
                cut_indices = np.where(energies < energy_cutoff)
                times = np.delete(times, cut_indices)
                energies = np.delete(energies, cut_indices)
            else:
                print('Potentially inaccurate large plastic calibration, beware radon washout!', file=datetime_logs)
            # Makes one bin for every ten seconds of the day
            if detector.processed:
                times = times - first_sec

            bins10sec = np.linspace(0.0, 86400.0, num=8641)

            # Creates numerical values for histograms using numpy
            hist_allday, bins_allday = np.histogram(times, bins=bins10sec)

            # Calculates mean and z-scores
            hist_allday_nz = hist_allday[hist_allday.nonzero()]
            mue = sum(hist_allday_nz) / len(hist_allday_nz)  # Mean number of counts/bin
            sigma = np.sqrt(mue)  # Standard deviation of the distribution

            # Removes outliers (3 sigma) and recalculates mue
            abs_zscores = np.abs((hist_allday_nz - mue) / sigma)
            hist_allday_nz = hist_allday_nz[np.argsort(abs_zscores)]
            hist_allday_nz = hist_allday_nz[::-1]
            zscores = np.sort(abs_zscores)[::-1]
            for k in range(len(abs_zscores)):
                if abs_zscores[k] < 3:
                    hist_allday_nz = hist_allday_nz[k:]
                    break

            mue = sum(hist_allday_nz) / len(hist_allday_nz)
            sigma = np.sqrt(mue)

            z_scores = np.array([])  # The z-scores themselves
            z_flags = np.array([])
            p = 0
            # Flags only those z-scores > 5
            for j in range(len(hist_allday)):
                if hist_allday[j] > mue:  # Peak
                    z_scores = np.append(z_scores, ((hist_allday[j] - mue) / sigma))
                    if z_scores[j] >= 5:
                        z_flags = np.append(z_flags, j)
                        p += 1
                elif hist_allday[j] < mue:  # Valley
                    z_scores = np.append(z_scores, ((hist_allday[j] - mue) / sigma))

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
                if glow_length > 0 and bins10sec[flag] - 10 == previous_time:
                    glow_length += 1
                if glow_length > 0 and bins10sec[flag] - 10 > previous_time:
                    glow = sc.PotentialGlow(glow_start, glow_length)
                    potential_glow_list.append(glow)
                    glow_start = flag
                    glow_length = 1
                if flag == z_flags[-1]:
                    glow = sc.PotentialGlow(glow_start, glow_length)
                    potential_glow_list.append(glow)

                previous_time = bins10sec[flag]

            sm.print_logger('Done.', datetime_logs)
            if len(potential_glow_list) == 0:
                sm.print_logger(f'There were no potential glows for the date {date_timestamp}', datetime_logs)
            else:
                # Logs potential glows and sorts them in descending order depending on their highest z-score
                highest_scores = []
                print('\n', file=datetime_logs)
                print('Potential glows:', file=datetime_logs)
                for glow in potential_glow_list:
                    highest_score = glow.highest_zscore(z_scores)
                    highest_scores.append(highest_score)
                    beginning, length = glow.glow_length_and_beginning_seconds(bins10sec)
                    print(f'{datetime.datetime.utcfromtimestamp(beginning + first_sec)} UTC ({beginning} '
                          f'seconds of day), {length} seconds long, highest z-score: {highest_score}',
                          file=datetime_logs)

                glow_sorting_order = np.argsort(highest_scores)
                glow_sorting_order = glow_sorting_order[::-1]
                potential_glow_list = [potential_glow_list[s] for s in glow_sorting_order]

                # Establishes detector location based on year and detector name
                location = sm.location(unit, YEAR)

                # Plotting the histograms
                sm.print_logger('\n', datetime_logs)
                sm.print_logger('Generating Histogram...', datetime_logs)
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
                n, bins, patches = ax1.hist(times, bins=bins10sec, alpha=0.5, color='r')
                ax1.set_xlabel('Seconds of Day (UT)')
                ax1.set_ylabel('Counts/10-second bin')
                ax1.axhline(y=(mue + 5 * sigma), color='blue', linestyle='dashed', linewidth=2)

                # Creates legend
                allday_data = mpatches.Patch(color='r', label='All Energies')
                allday_5sigma = mpatches.Patch(color='blue', label='5 Sigma Above All Energies', linestyle='dashed')
                ax1.legend(handles=[allday_data, allday_5sigma, ], bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
                ax1.grid(True)

                for a in range(4):
                    try:
                        glow = potential_glow_list[a]
                        sm.hist_subplotter(ax_list[a], glow, times, bins10sec, mue, sigma)
                    except IndexError:
                        continue

                plt.tight_layout()
                plt.plot()

                # Saves the histograms:
                hist_path = f'{sm.results_loc()}Results/{unit}/{full_day_string}/'
                sm.path_maker(hist_path)
                plt.savefig(f'{hist_path}{full_day_string}_histogram.png', dpi=500)
                plt.close(figu)
                sm.print_logger('Done', datetime_logs)
                sm.print_logger('\n', datetime_logs)

            datetime_logs.close()
