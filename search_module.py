"""Functions for use in the TGF and glow search program.

This module contains all the custom functions that are used by the TGF and glow search program. It is also meant to
serve as a centralized location for file paths, detector locations, and other relevant information.

Functions:
    C_raw_data_loc:
        Returns the custom data path specified by the user.
    results_loc:
        Returns the path where all the program's results (i.e. the log, scatter plots, etc.) are saved.
    G_raw_data_loc:
        Returns the path where GODOT raw data is stored on Sol.
    T_raw_data_loc:
        Returns the path where THOR raw data is stored on Sol.
    S_raw_data_loc:
        Returns the path where Santis instrument raw data is stored on Sol.
    CR_raw_data_loc:
        Returns the path where Croatia instrument raw data is stored on Sol.
    G_processed_data_loc:
        Returns the path where GODOT processed data is stored on Sol.
    memory_allowance:
        Returns maximum amount of memory that the program can use (in bytes), not including data.
    T_eRC:
        Returns the eRC serial number for the requested THOR unit.
    location:
        Returns the location of the requested instrument at the specified time.
    days_per_month:
        Returns the number of days in a given month.
    print_logger:
        Prints the given string in both the console and in the specified text file.
    chunk_unpickler:
        Unpickles and loads daily chunks for the program's low memory mode.
    path_maker:
        Makes the directory path specified by a string.
    channel_to_mev:
        Converts all the values of an array containing photon energies from channel voltage to MeV based on two
        conversion factors obtained from the associated scintillator's calibration.

"""

import os as os
import pickle as pickle


# Info functions:

def C_raw_data_loc():
    """Custom raw data location."""
    # Enter the custom raw data filepath inbetween the quotes.
    path = ''
    return path


def results_loc():
    """Location where results from search.py are exported (same as program location by default)."""
    path = ''
    return path


def G_raw_data_loc():
    """Location of GODOT raw data on Sol."""
    path = '/media/AllDetectorData/Detectors/GODOT/Data'
    return path


def T_raw_data_loc():
    """Location of THOR raw data on Sol."""
    path = '/media/AllDetectorData/Detectors/THOR'
    return path


def S_raw_data_loc():
    """Location of Santis instrument raw data on Sol."""
    path = '/media/AllDetectorData/Detectors/SANTIS/Data'
    return path


def CR_raw_data_loc():
    """Location of Croatia instrument raw data on Sol."""
    path = '/media/AllDetectorData/Detectors/SANTIS/Data'
    return path


def G_processed_data_loc():
    """Location of GODOT processed data on Sol."""
    path = '/media/godot/godot/monthly_processed'
    return path


def memory_allowance():
    """Returns maximum amount of memory that the program is allowed to use (in bytes), not including data."""
    memory = 100e6
    return memory


def T_eRC(unit, day):
    """Returns a list of all THOR eRC serial numbers for the requested THOR unit."""
    # All lists are in this form: NaI, small plastic, medium plastic, large plastic
    THOR1 = ['4179', '4194', '4189', '4195']
    THOR2 = ['4182', '4172', '4167', '4187']
    THOR3 = ['4169', '4175', '4174', '4185']
    THOR4 = ['4177', '4191', '4192', '4181']
    THOR5 = ['4188', '4190', '4169' if int(day) >= 221109 else '4178', '4173']  # THOR5 MP was replaced on Nov. 9th 2022
    THOR6 = ['4186', '4176', '4183', '4180']

    THOR_dict = {'THOR1': THOR1, 'THOR2': THOR2, 'THOR3': THOR3, 'THOR4': THOR4, 'THOR5': THOR5, 'THOR6': THOR6}
    return THOR_dict[unit]


def location(unit, YEAR):
    """Returns the location of a detector for the specified date."""
    # To add another detector/location simply add another if/elif statement
    # In the future it might be better to just import all of these from a text file or something
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
        location = 'at no location listed'

    return location


def days_per_month(month, year):
    """Returns the number of days in the requested month based on the year."""
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
def print_logger(string, logfile):
    """Prints the specified string to both the console and the specified text file."""
    print(string)
    print(string, file=logfile)


def chunk_unpickler(chunk_path):
    """Unpickles and loads daily chunks for the program's low memory mode."""
    chunk_file = open(chunk_path, 'rb')
    chunk = pickle.load(chunk_file)
    chunk_file.close()
    return chunk


def path_maker(path):
    """Checks to see if a directory path corresponding to the given string exists and, if not, creates it."""
    if not os.path.exists(path):
        os.makedirs(path)


def channel_to_mev(energy_array, channel, scintillator):
    """Uses compton edges/photo peaks obtained from scintillator calibration to convert energy channels into MeV

    Parameters
    ----------
    energy_array : np.array
        The array containing all the energies for either the large plastic or sodium iodide scintillator.
    channel : np.array
        An array containing the energy channels corresponding to the compton edges/photo peaks.
    scintillator : str
        The label corresponding to the scintillator whose energies are being converted.

    Returns
    -------
    np.array
        An array full of float values. These are the energies in MeV.

    """
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
