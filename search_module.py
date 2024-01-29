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
        Returns the calendar date after the one given as an argument (in yymmdd format).
    roll_date_backward:
        Returns the calendar date before the one given as an argument (in yymmdd format).
    full_date_to_short:
        Converts a date string of the form yyyy-mm-dd to the form yymmdd.
    short_to_full_date:
        Converts a date string of the form yymmdd to the form yyyy-mm-dd.
    print_logger:
        Prints the given string in both the console and in the specified text file.
    chunk_unpickler:
        Unpickles and loads daily chunks for the program's low memory mode.
    path_maker:
        Makes the directory path specified by a string.
    channel_to_mev:
        Converts all the values of an array containing photon energies from channel voltage to MeV based on two
        conversion factors obtained from the associated scintillator's calibration.
    convert_to_local:
        Converts the detector date and event time to what they would actually be in local time.
    get_weather_conditions:
        Scrapes weather underground and returns the weather at the approximate time of an event.
    scrape_weather:
        Scrapes weather from weather underground and returns the results as a pandas data frame.
    dst_status:
        Returns string statuses depending on whether a day falls inside/outside/on the edge of dst.
    dst_conversion:
        Returns an updated utc to local conversion number depending on the given date and time.
    convert_clock_hour:
        Converts a timestamp of the form hh:mm AM/PM into seconds since the beginning of the day.
    weather_from_code:
        Returns the weather for each code given by the function get_weather_conditions.

"""

import os as os
import contextlib as contextlib
import pickle as pickle
import datetime as dt
import glob as glob
import json as json
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
import pandas as pd

# Constants needed by various functions
TWO_AM = 7200  # Number of seconds of the day corresponding to 2:00AM
SEC_PER_HOUR = 3600
SEC_PER_DAY = 86400
CENTURY = '20'


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
    return '/media/AllDetectorData/Detectors/GODOT/Data'


def T_raw_data_loc():
    """Location of THOR raw data on Sol."""
    return '/media/AllDetectorData/Detectors/THOR'


def S_raw_data_loc():
    """Location of Santis instrument raw data on Sol."""
    return '/media/AllDetectorData/Detectors/SANTIS/Data'


def CR_raw_data_loc():
    """Location of Croatia instrument raw data on Sol."""
    return '/media/AllDetectorData/Detectors/SANTIS/Data'


def G_processed_data_loc():
    """Location of GODOT processed data on Sol."""
    return '/media/godot/godot/monthly_processed'


def memory_allowance():
    """Returns maximum amount of memory that the program is allowed to use (in bytes), not including data."""
    return 100e6


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


def location(unit, date):
    """Returns the location of a detector for the specified date."""
    if unit == 'GODOT':
        deployment_file_loc = G_raw_data_loc()[:-5]
    elif unit[0:4] == 'THOR':
        deployment_file_loc = T_raw_data_loc()
    else:
        deployment_file_loc = S_raw_data_loc()[:-5]

    for file in glob.glob(f'{deployment_file_loc}/{unit.lower()}_*_*.json'):
        if int(file[6:12]) <= date <= int(file[13:19]):
            with open(file, 'r') as deployment:
                return json.load(deployment)

    return {'Location': 'no location listed', 'Instrument': unit, 'Start date': '', 'End date': '',
            'UTC conversion to local time': '', 'Nearest weather station': '', 'Daylight Savings?': '',
            'Latitude (N)': '', 'Longitude (E, 0-360)': '', 'Altitude (km)': '',
            'Notes': ''}


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
    """Returns the calendar date after the one given as an argument. (in yymmdd format)"""
    date_int = int(date_str)
    date_int += 1
    date_str = str(date_int)
    # Month rollover
    if int(date_str[4:]) > days_per_month(int(date_str[2:4]), int(CENTURY + date_str[0:2])):
        date_int = date_int + 100 - (int(date_str[4:]) - 1)
        date_str = str(date_int)

    # Year rollover
    if int(date_str[2:4]) > 12:
        date_int = (date_int // 10000 + 1) * 10000 + 101
        date_str = str(date_int)

    return date_str


def roll_date_backward(date_str):
    """Returns the calendar date before the one given as an argument (in yymmdd format)."""
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
        date_int += days_per_month(int(str(date_int)[2:4]), int(CENTURY + date_str[0:2]))
        date_str = str(date_int)

    return date_str


def full_date_to_short(full_date_str):
    """Converts a date string of the form yyyy-mm-dd to the form yymmdd."""
    return full_date_str[2:].replace('-', '')


def short_to_full_date(date_str):
    """Converts a date string of the form yymmdd to the form yyyy-mm-dd."""
    return f'{CENTURY}{date_str[0:2]}-{date_str[2:4]}-{date_str[4:]}'


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


def convert_to_local(detector, event_time):
    """Converts the detector date and event time to what they would actually be in local time.

    Parameters
    ----------
    detector : sc.Detector
        The detector object where the date is stored
    event_time : float
        The time of the day when the event occurred in seconds since the beginning of the day.

    Returns
    -------
    str / float
        An updated date in local time and an updated event time in local time.

    """

    date_str = detector.date_str
    timezone_conversion = detector.location['UTC conversion to local time']

    # Just in case the event happened in the ~300 seconds of the next day typically included in the dataset
    if event_time > SEC_PER_DAY:
        event_time -= SEC_PER_DAY
        date_str = roll_date_forward(date_str)

    # Corrects the UTC conversion if we're in daylight savings time
    if detector.location['Daylight savings?'].lower() in ['yes', 'true']:  # Not sure which we're using yet
        timezone_conversion = dst_conversion(date_str, event_time, timezone_conversion)

    # If the event happened the next day local time
    if (event_time + (SEC_PER_HOUR * timezone_conversion)) > SEC_PER_DAY:
        date_str = roll_date_forward(date_str)
        event_time = (event_time + (SEC_PER_HOUR * timezone_conversion)) - SEC_PER_DAY
    # If the event happened the previous day local time
    elif (event_time + (SEC_PER_HOUR * timezone_conversion)) < 0:
        date_str = roll_date_backward(timezone_conversion)
        event_time = (event_time + (SEC_PER_HOUR * timezone_conversion)) + SEC_PER_DAY
    else:
        event_time = event_time + (SEC_PER_HOUR * timezone_conversion)

    return short_to_full_date(date_str), event_time


def get_weather_conditions(full_date_str, event_time, detector, weather_cache):
    """Scrapes weather underground and returns the weather at the approximate time of an event.

    Parameters
    ----------
    full_date_str : str
        The date that the event occurred on (in local time) in yyyy-mm-dd format.
    event_time : float
        The time that the event occurred at during the day (in local time) in units of seconds since beginning of day.
    detector : sc.Detector
        The detector object that contains the name of the nearest weather station.
    weather_cache : dict
        A cache containing weather tables that have already been retrieved. Keys are dates in yyyy-mm-dd format.

    Returns
    -------
    int
        A score corresponding to the weather conditions around the time of the event. See the function
        weather_from_score for a summary of what the scores mean.

    """
    if full_date_str in weather_cache and weather_cache[full_date_str] is not None:
        weather_table = weather_cache[full_date_str]
    else:
        weather_table = scrape_weather(full_date_str, detector.location['Station'])
        weather_cache[full_date_str] = weather_table

    if weather_table is not None:
        # Finds the time in the table that's closest to the time of the event
        index = 0
        best_diff = float('inf')
        best_index = 0
        for clock_hour in weather_table['Time']:
            if type(clock_hour) != float:
                time_sec = convert_clock_hour(clock_hour)
                time_diff = abs(event_time - time_sec)
                if time_diff < best_diff:
                    best_diff = time_diff
                    best_index = index
            else:
                break

            index += 1

        # Gets the weather conditions at the closest hour to the event and the surrounding hour_padding hours
        weather = []
        hour_padding = 3
        for i in range(best_index - hour_padding, best_index + hour_padding + 1):
            if 0 <= i < index:
                weather.append(weather_table['Condition'][i])
            else:
                weather.append(None)

        heavy_rain = False
        rain = False
        for condition in weather:
            if condition:
                for variation in ['Thunder', 'T-Storm', 'Storm', 'Lightning', 'Hail']:
                    if variation in condition:
                        return 1

                if 'Heavy' in condition:
                    heavy_rain = True
                elif 'Rain' in condition:
                    rain = True

        if heavy_rain:
            return 0.75
        elif rain:
            return 0.5

        return 0
    else:
        return -1


def scrape_weather(full_date_str, station):
    """Scrapes weather from weather underground and returns the results as a pandas data frame.

    Parameters
    ----------
    full_date_str : str
        The date that weather data is being requested for in yyyy-mm-dd format.
    station : str
        The four-letter name of the weather station that data is being requested for.

    Returns
    -------
    pd.Dataframe
        A pandas dataframe with weather information for the specified day.

    """
    try:
        # Note: selenium and lxml modules are required to make this work. Install them
        chrome_options = Options()
        chrome_options.add_argument('--headless=new')  # Runs the chrome client in headless mode
        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):  # Prevents selenium from printing status stuff
            driver = webdriver.Chrome(options=chrome_options)

            url = f'https://www.wunderground.com/history/daily/{station}/date/{full_date_str}'

            driver.get(url)
            tables = WebDriverWait(driver, 20).until(ec.presence_of_all_elements_located((By.CSS_SELECTOR, "table")))

        table = pd.read_html(tables[1].get_attribute('outerHTML'))[0]

        return table  # This is a dataframe containing the table we want

    except:
        return None


def dst_status(date_str):
    """Returns string statuses depending on whether a day falls inside/outside/on the edge of dst.

    Parameters
    ----------
    date_str : str
        The date to be checked in yymmdd format.

    Returns
    -------
    str
        A status for the date: inside if the date is inside dst, outside if out, or beginning/end for the boundaries.

    """
    year = int(CENTURY + date_str[0:2])
    month = int(date_str[2:4])
    day = int(date_str[4:])

    # January, February, and December are never DST
    if month < 3 or month > 11:
        return 'outside'
    # April to October are always DST
    elif 3 < month < 11:
        return 'inside'
    # DST starts on the second Sunday of March (which is always between the 8th and the 14th)
    elif month == 3:
        second_sunday = 8 + (6 - dt.datetime(year, month, 8).weekday())
        if day < second_sunday:
            return 'outside'
        elif day > second_sunday:
            return 'inside'
        else:
            return 'beginning'
    # DST ends on the first Sunday of November (so the previous Sunday must be before the 1st)
    else:
        first_sunday = 1 + (6 - dt.datetime(year, month, 1).weekday())
        if day < first_sunday:
            return 'inside'
        elif day > first_sunday:
            return 'outside'
        else:
            return 'end'


def dst_conversion(date_str, event_time, timezone_conversion):
    """Returns an updated utc to local conversion number depending on the given date and time.

    Parameters
    ----------
    date_str : str
        The date to be converted in yymmdd format.
    event_time : float
        The time that the event occurred in units of seconds since the beginning of the day.
    timezone_conversion : int
        A number giving the hour difference between local time and UTC.

    Returns
    -------
    int
        An updated timezone conversion that accounts for dst.

    """

    temp_time = event_time + (timezone_conversion * SEC_PER_HOUR)
    if temp_time > SEC_PER_DAY:
        temp_time -= SEC_PER_DAY
        temp_date = roll_date_forward(date_str)
    elif temp_time < 0:
        temp_time += SEC_PER_DAY
        temp_date = roll_date_backward(date_str)
    else:
        temp_date = date_str

    temp_date_status = dst_status(temp_date)
    if temp_date_status == 'inside':  # Squarely inside dst
        return timezone_conversion + 1
    elif temp_date_status == 'outside':  # Squarely outside dst
        return timezone_conversion
    elif temp_date_status == 'beginning':  # Beginning of dst (2nd Sunday of March at 2:00AM)
        if temp_time >= TWO_AM:
            return timezone_conversion + 1
        else:
            return timezone_conversion
    else:  # End of dst (1st Sunday of November at 2:00AM)
        if (temp_time + SEC_PER_HOUR) >= TWO_AM:  # + sec_per_hour b/c temp time should be in dst
            return timezone_conversion
        else:
            return timezone_conversion + 1


def convert_clock_hour(clock_hour):
    """Converts a timestamp of the form hh:mm AM/PM into seconds since the beginning of the day."""

    meridiem = clock_hour.split()[1]
    hour = int(clock_hour.split()[0].split(':')[0])
    minute = int(clock_hour.split()[0].split(':')[1])

    # Converting from 12 hour time to 24 hour time
    if meridiem == 'AM' and hour == 12:  # midnight
        hour = 0
    elif meridiem == 'PM' and hour == 12:  # noon
        pass
    elif meridiem == 'PM':  # PM conversion
        hour += 12

    return float((hour * 3600) + (minute * 60))


def weather_from_score(score):
    """Returns the weather for each code given by the function get_weather_conditions."""
    if score == 0:
        return 'fair'
    elif score == 0.5:
        return 'light rain'
    elif score == 0.75:
        return 'heavy rain'
    elif score == 1:
        return 'Lightning or hail'
    else:
        return 'error getting weather data'
