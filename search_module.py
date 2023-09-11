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
    roll_date_forward:
        Returns the calendar date after the one given as an argument.
    roll_date_backward:
        Returns the calendar date before the one given as an argument.
    print_logger:
        Prints the given string in both the console and in the specified text file.
    chunk_unpickler:
        Unpickles and loads daily chunks for the program's low memory mode.
    path_maker:
        Makes the directory path specified by a string.
    channel_to_mev:
        Converts all the values of an array containing photon energies from channel voltage to MeV based on two
        conversion factors obtained from the associated scintillator's calibration.
    utc_to_local_date:
        Converts the detector date to the correct local one based on the time of the supplied event.
    scrape_weather:
        Scrapes weather from weather underground and returns the results as a pandas data frame.
    convert_timestamp:
        Converts a timestamp of the form hh:mm AM/PM into seconds since the beginning of the day.
    get_weather:
        Scrapes weather underground and returns the weather at the approximate time of an event.

"""

import os as os
import contextlib as contextlib
import pickle as pickle
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd


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


def roll_date_forward(date_str):
    """Returns the calendar date after the one given as an argument."""
    century = '20'
    date_int = int(date_str)
    date_int += 1
    date_str = str(date_int)
    # Month rollover
    if int(date_str[4:]) > days_per_month(int(date_str[2:4]), int(century + date_str[0:2])):
        date_int = date_int + 100 - (int(date_str[4:]) - 1)
        date_str = str(date_int)

    # Year rollover
    if int(date_str[2:4]) > 12:
        date_int = (date_int // 10000 + 1) * 10000 + 101
        date_str = str(date_int)

    return date_str


def roll_date_backward(date_str):
    """Returns the calendar date before the one given as an argument."""
    century = '20'
    date_int = int(date_str)
    date_int -= 1
    date_str = str(date_int)
    # Year rollback
    if int(date_str[2:]) == 100:  # This would be January 0th because int(0100) = 100
        date_int = (date_int // 10000 - 1) * 10000 + 1231  # December 31st of the previous year
        date_str = str(date_int)

    # Month rollback
    if int(date_str[4:]) == 0:
        date_int -= 100
        date_int += days_per_month(int(str(date_int)[2:4]), int(century + date_str[0:2]))
        date_str = str(date_int)

    return date_str


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


def utc_to_local_date(detector, event_time):
    """Converts the detector date to the correct local one based on the time of the supplied event."""
    utc_conversion = int(detector.location['UTC_conversion'])

    sec_per_hour = 3600
    sec_per_day = 86400
    century = '20'

    if (event_time + (sec_per_hour * utc_conversion)) > sec_per_day:  # If the event happened the next day local time
        full_day_string = roll_date_forward(detector.full_day_string)
    elif (event_time + (sec_per_hour * utc_conversion)) < 0:
        full_day_string = roll_date_backward(detector.full_day_string)  # If the event happened the previous day lt
    else:
        return detector.date_timestamp

    return f'{century}{full_day_string[0:2]}-{full_day_string[2:4]}-{full_day_string[4:]}'


def scrape_weather(date, station):
    """Scrapes weather from weather underground and returns the results as a pandas data frame.

    Parameters
    ----------
    date : str
        The date to scrape weather for in the format yyyy-mm-dd.
    station : str
        The weather station to scrape weather info from.

    Returns
    -------
    pd.DataFrame
        A pandas data frame containing the weather info for the requested date.

    """
    try:
        # Note: selenium and lxml modules are required to make this work. Install them
        chrome_options = Options()
        chrome_options.add_argument('--headless=new')  # Runs the chrome client in headless mode
        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):  # Prevents selenium from printing status stuff
            driver = webdriver.Chrome(options=chrome_options)

            url = f'https://www.wunderground.com/history/daily/{station}/date/{date}'

            driver.get(url)
            tables = WebDriverWait(driver, 20).until(EC.presence_of_all_elements_located((By.CSS_SELECTOR, "table")))

        return pd.read_html(tables[1].get_attribute('outerHTML')[0])  # This is a dataframe containing the table we want

    except:
        return None


def convert_timestamp(timestamp):
    """Converts a timestamp of the form hh:mm AM/PM into seconds since the beginning of the day.

    Parameters
    ----------
    timestamp : str
        A string containing the timestamp to be converted (in the format hh:mm AM/PM).

    Returns
    -------
    float
        The number of seconds of the day corresponding to the timestamp.

    """
    meridiem = timestamp.split()[1]
    hour = int(timestamp.split()[0].split(':')[0])
    minute = int(timestamp.split()[0].split(':')[1])

    # Converting from 12 hour time to 24 hour time
    if meridiem == 'AM' and hour == 12:  # midnight
        hour = 0
    elif meridiem == 'PM' and hour == 12:  # noon
        pass
    elif meridiem == 'PM':  # PM conversion
        hour += 12

    return float((hour * 3600) + (minute * 60))


def get_weather(detector, times, event):
    """Scrapes weather underground and returns the weather at the approximate time of an event.

    Parameters
    ----------
    detector : Detector
        The detector object being used in the search.
    times : np.array
        A numpy array containing the list mode times that the event was found in.
    event : ShortEvent
        A short event object corresponding to the event.

    Returns
    -------
    str
        A string containing the weather conditions at the approximate time of the event.

    """
    event_time = times[event.start] - 86400 if times[event.start] > 86400 else times[event.start]

    date = utc_to_local_date(detector, event_time)
    station = detector.location['nearest_station']
    weather_table = scrape_weather(date, station)

    if weather_table is not None:
        # Finds the time in the table that's closest to the time of the event
        index = 0
        best_diff = float('inf')
        best_index = 0
        for time in weather_table['Time']:
            if type(time) != float:
                time_sec = convert_timestamp(time)
                time_diff = abs(event_time - time_sec)
                if time_diff < best_diff:
                    best_diff = time_diff
                    best_index = index
                else:
                    break

            index += 1

        # Returns the weather condition at the closest time and the surrounding two hours
        weather_string = ''
        if (best_index - 1) >= 0:
            weather_string += f"{weather_table['Condition'][best_index - 1]} (hr before), "

        weather_string += f"{weather_table['Condition'][best_index]} (during)"

        if best_index < index:
            weather_string += f", {weather_table['Condition'][best_index + 1]} (hr after)"

        return weather_string

    else:
        return 'not found'
