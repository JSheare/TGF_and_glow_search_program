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
def channel_to_mev(energy_array, channel, scintillator):
    if scintillator == 'NaI':
        K40 = 1.46  # Photo-peak photon energy for Potassium 40 (MeV)
        T = 2.60  # Photo-peak photon energy for Thorium (MeV)
    else:  # i.e. plastic scintillators
        K40 = 1.242  # Compton edge photon energy for Potassium 40 (MeV)
        T = 2.381  # Compton edge photon energy for Thorium (MeV)

    channel1 = channel[0]
    channel2 = channel[1]
    a = (K40 - T)/(channel1-channel2)
    b = T - a*channel2
    energy_array = a*energy_array + b
    return energy_array


# Makes the flagged z-score subplots for the glow search histograms
def hist_subplotter(ax, glow, times, bins10sec, mue, sigma):
    padding = 20
    c = 0 if (glow.peak_index - padding) < 0 else (glow.peak_index - padding)
    d = 8670 if (glow.peak_index + padding) > 8670 else (glow.peak_index + padding)

    subbins = bins10sec[c:d]
    ax.hist(times, bins=subbins, alpha=0.5, color='c')
    ax.set_xlabel('Seconds of Day (UT)')
    ax.set_ylabel('Counts/10-sec', )
    ax.axhline(y=(mue + 5 * sigma), color='blue', linestyle='dashed', linewidth=2)
    ax.grid(True)
