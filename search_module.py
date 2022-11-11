import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Functions used in search.py. Also, a centralized location for file paths and detector locations


# Info functions:

# Custom raw data location
def C_raw_data_loc():
    # Enter the custom raw data filepath inbetween the quotes.
    path = ''
    return path


# Location where results from search.py are exported (same as program location by default)
def results_loc():
    path = ''
    return path


# Location of GODOT raw data
def G_raw_data_loc():
    path = '/media/AllDetectorData/Detectors/GODOT/Data'
    return path


# Location of THOR raw data
def T_raw_data_loc():
    path = '/media/AllDetectorData/Detectors/THOR'
    return path


# Location of Santis instrument raw data
def S_raw_data_loc():
    path = '/media/AllDetectorData/Detectors/SANTIS/Data'
    return path


# Location of GODOT processed data
def G_processed_data_loc():
    path = '/media/godot/godot/monthly_processed'
    return path


# THOR eRC values:
def T_eRC(unit):
    # all lists are in this form: NaI, small plastic, medium plastic, large plastic
    THOR1 = ['4179', '4194', '4189', '4195']
    THOR2 = ['4182', '4172', '4167', '4187']
    THOR3 = ['4169', '4175', '4174', '4185']
    THOR4 = ['4177', '4191', '4192', '4181']
    THOR5 = ['4188', '4190', '4178', '4173']
    THOR6 = ['4186', '4176', '4183', '4180']

    THOR_dict = {'THOR1': THOR1, 'THOR2': THOR2, 'THOR3': THOR3, 'THOR4': THOR4, 'THOR5': THOR5, 'THOR6': THOR6}
    return THOR_dict[unit]


# Detector Locations. To add another detector/location simply add another if/elif statement
def location(unit, YEAR):
    if unit == 'GODOT':
        if YEAR == 2015:
            location = 'at HAWC, Mexico until September, then at Uchinada'
        elif YEAR == 2016:
            location = 'at HAWC, Mexico'
        elif YEAR == 2017:
            location = 'at HAWC, Mexico'
        elif YEAR == 2018:
            location = 'at Uchinada, Japan'
        elif YEAR == 2019:
            location = 'at Uchinada, Japan'
        elif YEAR == 2020:
            location = 'at Uchinada, Japan'
        else:  # YEAR == 2021-20xx
            location = 'at Uchinada, Japan'

    elif unit[0:4] == 'THOR':  # Temporary. Eventually this will contain a list of all THOR unit locations as well
        if unit == 'THOR5':
            location = 'at mount Fuji, Japan'  # temp
        else:
            location = 'at no listed location'

    elif unit == 'SANTIS':
        location = 'Santis, Switzerland'
    else:
        location = 'at no listed location'

    return location


# Returns the number of days in the requested month
def days_per_month(month, year):
    month_dict = {'1': 31,  # January
                  '2': 29 if year % 4 == 0 or (year % 100 != 0 and year % 400 == 0) else 28,  # February
                  '3': 31,  # March
                  '4': 30,  # April
                  '5': 31,  # May
                  '6': 30,  # June
                  '7': 31,  # July
                  '8': 31,  # August
                  '9': 30,  # September
                  '10': 31,  # October
                  '11': 30,  # November
                  '12': 31}  # December
    return month_dict[str(month)]


# Program functions:
# Prints the specified string to both the console and the log
def print_logger(string, logfile):
    print(string)
    print(string, file=logfile)


# Checks to see if a given path exists and, if not, creates it
def path_maker(path):
    if not os.path.exists(path):
        os.makedirs(path)


# Uses compton edges obtained from scintillator calibration to convert energy array values into actual energies in MeV
def channel_to_mev(energy_array, channel1, channel2, eRC):  # review this
    if eRC == '1491':
        K40 = 9.982624235104672506e-13  # Compton edge photon wavelength for Potassium 40 (meters)
        T = 5.207231961360773518e-13  # Compton edge photon wavelength for Thorium (meters)
    else:
        K40 = 8.492068013698631405e-13  # Photo-peak photon wavelength for Potassium 40 (meters)
        T = 4.768622807692308078e-13  # Photo-peak photon wavelength for Thorium (meters)
    h = 4.1357e-21  # Planck's constant in MeV*s

    a = (K40 - T) / (channel1 - channel2)
    b = K40 - (channel1 * (K40 - T)) / (channel1 - channel2)
    energy_array = (h * 2.998e8) / (a * energy_array + b)
    return energy_array


# Makes the flagged z-score subplots for the glow search histograms
def hist_subplotter(ax, glow, times, bins10sec, mue, sigma):
    padding = 20
    c = 0 if (glow.peak_index - padding) < 0 else (glow.peak_index - padding)
    d = 8639 if (glow.peak_index + padding) > 8639 else (glow.peak_index + padding)

    subbins = bins10sec[c:d]
    ax.hist(times, bins=subbins, alpha=0.5, color='c')
    ax.set_xlabel('Seconds of Day (UT)')
    ax.set_ylabel('Counts/10-sec', )
    ax.axhline(y=(mue + 5 * sigma), color='blue', linestyle='dashed', linewidth=2)
    ax.grid(True)


# Makes the energy spectra histograms for each GODOT scintillator
def g_spectra_maker(energies, date_timestamp, full_day_string, detector, eRC):  # review this
    plt.figure(figsize=[20, 11.0])
    energy_bins = np.linspace(0.0, 15008.0, num=939)
    energy_hist, bin_edges = np.histogram(energies, bins=energy_bins)
    flagged_indices = np.array([])

    if eRC == '1491':
        scintillator = 'Large Plastic'
        window_size = 19
        poly_order = 2
        curve_change = 1000
        # Makes a rough "response function" for plotting purposes.
        smoothed_energies = signal.savgol_filter(energy_hist, window_size, poly_order)  # Savitzky-Golay filter
        plt.plot(energy_bins[0:938], smoothed_energies, color='green', alpha=0.75)

        # Finds all the points on the response function where the concavity changes from positive to negative
        second_derivative = signal.savgol_filter(energy_hist, window_size, poly_order, 2)
        curvature_sign = np.sign(second_derivative)
        for jk in range(len(second_derivative)):
            if jk > 0 and curvature_sign[jk] == -1 and curvature_sign[jk - 1] in [0, 1]:
                lower_limit = energy_hist[jk - int(((window_size - 1) / 2))]
                upper_limit = energy_hist[jk + int(((window_size - 1) / 2))]
                # Only flags those points whose surrounding curvature is sufficiently hilly
                if np.abs(lower_limit - upper_limit) > curve_change:
                    flagged_indices = np.append(flagged_indices, jk)
        if len(flagged_indices) < 2:
            raise Exception('One or more compton edges not found')

        if len(flagged_indices) > 2:
            raise Exception('Too many concavity changes found')

    else:
        scintillator = 'Sodium Iodide'
        # Takes the sum of each bin with its two closest neighboring bins on either side
        sums = energy_hist
        for i in range(2):
            sums += (np.roll(energy_hist, i+1) + np.roll(energy_hist, i-1))

        # Looks for the location of the maximum sum within the two bands where the peaks are likely to be
        band_starts = [38, 94]
        band_ends = [75, 125]
        for th in range(len(band_starts)):
            band_max = np.argmax(sums[band_starts[th]:band_ends[th]]) + int(band_starts[th])
            flagged_indices = np.append(flagged_indices, band_max)

    # Plots the actual spectrum
    plt.title(f'Energy Spectrum for {scintillator}, {str(date_timestamp)}', loc='center')
    plt.xlabel('Energy Channel')
    plt.ylabel('Counts/bin')
    plt.yscale('log')
    plt.hist(energies, bins=energy_bins, color='r', rwidth=0.5, zorder=1)

    # Saves energy bins corresponding to the desired energies and plots them as vertical lines
    energy1 = energy_bins[int(flagged_indices[0])]
    energy2 = energy_bins[int(flagged_indices[1])]
    plt.vlines(energy_bins[flagged_indices.astype(int)], 0, 1e6, zorder=2, alpha=0.75)

    # Saves the figure
    sp_path = f'{results_loc()}Results/{detector}/{full_day_string}/'
    path_maker(sp_path)
    plt.savefig(f'{sp_path}{scintillator}_Spectrum.png', dpi=500)
    plt.clf()

    return energy1, energy2
