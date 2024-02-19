"""Tools for use by the TGF search program and its modules."""
import os as os
import contextlib as contextlib
import pickle as pickle
import datetime as dt
import numpy as np
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec

import tgf_search.parameters as params


def print_logger(string, logfile):
    """Prints the specified string to both the console and the specified text file.

    Parameters
    ----------
    string : str
        The string to be printed/logged
    logfile : file
        The file where the string should be written.

    """

    print(string)
    print(string, file=logfile)


def make_path(path):
    """Checks to see if a directory path corresponding to the given string exists and, if not, creates it.

    Parameters
    ----------
    path : str
        The path to be created.

    """

    if not os.path.exists(path):
        os.makedirs(path)


def days_per_month(month, year):
    """Returns the number of days in the requested month based on the year.

    Parameters
    ----------
    month : int
        The month of the year (1-12).
    year : int
        The desired year to check.

    Returns
    -------
    int
        The number of days in the specified month for the specified year.

    """

    if month == 1:  # January
        return 31
    elif month == 2:  # February
        return 29 if year % 4 == 0 or (year % 100 != 0 and year % 400 == 0) else 28
    elif month == 3:  # March
        return 31
    elif month == 4:  # April
        return 30
    elif month == 5:  # May
        return 31
    elif month == 6:  # June
        return 30
    elif month == 7:  # July
        return 31
    elif month == 8:  # August
        return 31
    elif month == 9:  # September
        return 30
    elif month == 10:  # October
        return 31
    elif month == 11:  # November
        return 30
    else:  # December
        return 31


def roll_date_forward(date_str):
    """Returns the calendar date after the one given as an argument.

    Parameters
    ----------
    date_str : str
        The date to be rolled forward from (in yymmdd format).

    Returns
    -------
    str
        The calendar date after the one supplied (in yymmdd format).

    """

    date_int = int(date_str)
    date_int += 1
    date_str = str(date_int)
    # Month rollover
    if int(date_str[4:]) > days_per_month(int(date_str[2:4]), int(params.CENTURY + date_str[0:2])):
        date_int = date_int + 100 - (int(date_str[4:]) - 1)
        date_str = str(date_int)

    # Year rollover
    if int(date_str[2:4]) > 12:
        date_int = (date_int // 10000 + 1) * 10000 + 101
        date_str = str(date_int)

    return date_str


def roll_date_backward(date_str):
    """Returns the calendar date before the one given as an argument.

    Parameters
    ----------
    date_str : str
        The date to be rolled backward from (in yymmdd format).

    Returns
    -------
    str
        The calendar date before the one supplied (in yymmdd format).

    """
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
        date_int += days_per_month(int(str(date_int)[2:4]), int(params.CENTURY + date_str[0:2]))
        date_str = str(date_int)

    return date_str


def make_date_list(first_date, second_date):
    """Makes a list of dates from first_date to second_date (inclusive).

    Parameters
    ----------
    first_date : str
        The first date in the range.
    second_date : str
        The second date in the range.

    Returns
    -------
    list
        A list of dates on the specified range.

    """

    requested_dates = [first_date]
    if first_date != second_date:
        date_str = first_date
        while True:
            date_str = roll_date_forward(date_str)
            requested_dates.append(date_str)
            if date_str == second_date:
                break

    return requested_dates


def full_date_to_short(full_date_str):
    """Converts a date string of the form yyyy-mm-dd to the form yymmdd."""
    return full_date_str[2:].replace('-', '')


def short_to_full_date(date_str):
    """Converts a date string of the form yymmdd to the form yyyy-mm-dd."""
    return f'{params.CENTURY}{date_str[0:2]}-{date_str[2:4]}-{date_str[4:]}'


def get_first_sec(date_str):
    """Converts the given date string (in yymmdd format) to its first second in EPOCH time."""
    day = int(date_str[4:])
    month = int(date_str[2:4])
    year = int(params.CENTURY + date_str[0:2])
    return (dt.datetime(year, month, day, 0, 0) - dt.datetime(1970, 1, 1)).total_seconds()


def pickle_detector(detector, path_form):
    """Pickles detector objects.

    Parameters
    ----------
    detector : Detector object
        The detector to be pickled.
    path_form : str
        The name of the pickle file (including the directory path where it should be saved).

    """

    detector.log = None  # serializing open file objects results in errors
    detector.file_form = None  # serializing anonymous functions results in errors too
    with open(path_form, 'wb') as file:
        pickle.dump(detector, file)


def unpickle_detector(pickle_path, modes=None):
    """Unpickles detector objects.

    Parameters
    ----------
    pickle_path : str
        The path to the pickle file that the detector is stored in.
    modes : list
        Optional. A list of modes that the detector should operate under.

    Returns
    -------
    Detector
        A detector-type object.

    """

    if modes is None:
        modes = list()

    with open(pickle_path, 'rb') as file:
        detector = pickle.load(file)

        # Modes might not necessarily be the same in the serialized object
        detector.modes = modes
        detector.template = True if 'template' in modes else False
        detector.check_processed()
        detector.check_gui()
        return detector


def pickle_chunk(chunk, filename):
    """Pickles daily chunks for the program's low memory mode."""
    pickle_detector(chunk, filename)


def unpickle_chunk(chunk_path):
    """Unpickles and loads daily chunks for the program's low memory mode."""
    with open(chunk_path, 'rb') as file:
        chunk = pickle.load(file)

    return chunk


def filter_files(complete_filelist):
    """Returns an ordered list of files with duplicate/incompatible files filtered out."""
    unique_files = set()
    files = []
    extensions = []
    for file in complete_filelist:
        # Filters out trace mode files and .txtp files (whatever those are)
        if file[-3:] == 'xtr' or file[-4:] == 'txtp' or file[-5:] == 'xtrpp':
            continue
        # Filters out duplicates
        else:
            if file[-7:] == '.txt.gz':
                file = file[:-7]
                extension = '.txt.gz'
            elif file[-4:] == '.txt':
                file = file[:-4]
                extension = '.txt'
            else:
                file = file[:-4]
                extension = '.csv'

            if file not in unique_files:
                unique_files.add(file)
                files.append(file)
                extensions.append(extension)

    return [files[s] + extensions[s] for s in np.argsort(files)]  # Puts the files back in order


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
    if event_time > params.SEC_PER_DAY:
        event_time -= params.SEC_PER_DAY
        date_str = roll_date_forward(date_str)

    # Corrects the UTC conversion if we're in daylight savings time
    if detector.location['Daylight savings?'].lower() in ['yes', 'true']:  # Not sure which we're using yet
        timezone_conversion = dst_conversion(date_str, event_time, timezone_conversion)

    # If the event happened the next day local time
    if (event_time + (params.SEC_PER_HOUR * timezone_conversion)) > params.SEC_PER_DAY:
        date_str = roll_date_forward(date_str)
        event_time = (event_time + (params.SEC_PER_HOUR * timezone_conversion)) - params.SEC_PER_DAY
    # If the event happened the previous day local time
    elif (event_time + (params.SEC_PER_HOUR * timezone_conversion)) < 0:
        date_str = roll_date_backward(timezone_conversion)
        event_time = (event_time + (params.SEC_PER_HOUR * timezone_conversion)) + params.SEC_PER_DAY
    else:
        event_time = event_time + (params.SEC_PER_HOUR * timezone_conversion)

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
    year = int(params.CENTURY + date_str[0:2])
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
    """Returns an updated UTC to local conversion number depending on the given date and time.

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

    temp_time = event_time + (timezone_conversion * params.SEC_PER_HOUR)
    if temp_time > params.SEC_PER_DAY:
        temp_time -= params.SEC_PER_DAY
        temp_date = roll_date_forward(date_str)
    elif temp_time < 0:
        temp_time += params.SEC_PER_DAY
        temp_date = roll_date_backward(date_str)
    else:
        temp_date = date_str

    temp_date_status = dst_status(temp_date)
    if temp_date_status == 'inside':  # Squarely inside dst
        return timezone_conversion + 1
    elif temp_date_status == 'outside':  # Squarely outside dst
        return timezone_conversion
    elif temp_date_status == 'beginning':  # Beginning of dst (2nd Sunday of March at 2:00AM)
        if temp_time >= params.TWO_AM:
            return timezone_conversion + 1
        else:
            return timezone_conversion
    else:  # End of dst (1st Sunday of November at 2:00AM)
        if (temp_time + params.SEC_PER_HOUR) >= params.TWO_AM:  # + sec_per_hour b/c temp time should be in dst
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


def combine_data(detector):
    """Combines data from all scintillators into one set of arrays.

    Parameters
    ----------
    detector : Detector
        The detector object that data will be combined for

    Returns
    -------
        np.array
            times
                A numpy array containing a combined list of second-of-day times.
            energies
                A numpy array containing a combined list of energies.
            wallclock
                A numpy array containing a combined list of wallclock times.
            count_scints
                A numpy array containing a list of strings. Each entry corresponds to the scintillator
                that its corresponding count originated from.

    """

    times = []
    energies = []
    wallclock = []
    count_scints = []
    for scintillator in detector:
        times.append(detector.get_attribute(scintillator, 'time'))
        energies.append(detector.get_attribute(scintillator, 'energy'))
        wallclock.append(detector.get_attribute(scintillator, 'wc'))
        count_scints.append([scintillator] * len(times[-1]))

    times = np.concatenate(times)
    energies = np.concatenate(energies)
    wallclock = np.concatenate(wallclock)
    count_scints = np.concatenate(count_scints)

    sorting_order = np.argsort(times)
    times = times[sorting_order]
    energies = energies[sorting_order]
    wallclock = wallclock[sorting_order]
    count_scints = count_scints[sorting_order]
    return times, energies, wallclock, count_scints


def separate_data(times, energies, count_scints, start, stop):
    """Separates the combined data produced in combo mode into separate data for each scintillator, organized
    in two dictionaries: time_dict and energy_dict.

    Parameters
    ----------
    times : np.array
        A numpy array containing a combined list of second-of-day times.
    energies : np.array
        A numpy array containing a combined list of energies.
    count_scints : np.array
        A numpy array containing a list of strings. Each entry corresponds to the scintillator
        that its corresponding count originated from.
    start : int
        The beginning of the range to separate.
    stop : int
        The end of the range to separate.

    Returns
    -------
    dict
        time_dict
            A dictionary containing lists of count times for each scintillator on the specified range.
        energy_dict
            A dictionary containing lists of count energies for each scintillator on the specified range.

    """

    time_dict = dict()
    energy_dict = dict()
    for i in range(start, stop):
        scintillator = count_scints[i]
        if scintillator not in time_dict:
            time_dict[scintillator] = [times[i]]
            energy_dict[scintillator] = [energies[i]]
        else:
            time_dict[scintillator].append(times[i])
            energy_dict[scintillator].append(energies[i])

    # Converting back to numpy arrays
    for scintillator in time_dict:
        times = time_dict[scintillator]
        energies = energy_dict[scintillator]
        time_dict[scintillator] = np.array(times)
        energy_dict[scintillator] = np.array(energies)

    return time_dict, energy_dict


def channel_to_mev(energy_array, channels, scintillator):
    """Uses compton edges/photo peaks obtained from scintillator calibration to convert energy channels into MeV

    Parameters
    ----------
    energy_array : np.array
        The array containing all the energies for either the large plastic or sodium iodide scintillator.
    channels : np.array
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

    channel1 = channels[0]
    channel2 = channels[1]
    a = (K40 - T)/(channel1-channel2)
    b = T - a*channel2
    energy_array = a*energy_array + b
    return energy_array
