"""A script that automatically runs searches on unsearched data."""
import os as os
import subprocess as subprocess
import json as json
import glob as glob
import psutil as psutil
import platform as platform

import tgfsearch.tools as tl


def search(unit, checked_dates, data_path, results_path, auto_search_path):
    if unit not in checked_dates:
        checked_dates[unit] = []

    queue = []
    days = glob.glob(f'{data_path}/*')
    days.sort()
    for day in days:
        if day[-6:].isnumeric():
            has_data = True if len(os.listdir(day)) > 0 else False
            if has_data and day[-6:] not in checked_dates[unit]:
                queue.append(day[-6:])

    # Check the most recent stuff first
    queue = queue[::-1]

    # Note: script_path will probably need to be changed if the file structure of the package is ever modified
    script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '\\search.py'
    executable = 'python' if platform.system() == 'Windows' else 'python3'
    for day in queue:
        if unit == 'THOR1':
            # Control flow for THOR1 NOAA flights
            if 220930 <= int(day) <= 230208:
                command = [executable, script_path, day, day, 'THOR1', 'aircraft', 'skcali']
            else:
                command = [executable, script_path, day, day, 'THOR1']

        elif unit == 'SANTIS':
            # Control flow for SANTIS instrument becoming the CROATIA instrument
            if int(day) > 211202:
                command = [executable, script_path, day, day, 'CROATIA']
            else:
                command = [executable, script_path, day, day, 'SANTIS']

        else:
            command = [executable, script_path, day, day, f'{unit}']

        # Runs search.py for the day
        old_pwd = os.getcwd()
        os.chdir(results_path)
        result = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Uncomment this and comment the above line to get search.py console output (for console debugging)
        # result = subprocess.Popen(command)
        os.chdir(old_pwd)

        if result.wait() != 0:
            # Writes any errors that search.py runs into to text files so that they can be examined/fixed later
            tl.make_path(f'{auto_search_path}/Error Logs')
            with open(f'{auto_search_path}/Error Logs/{unit}_{day}_Error.txt', 'w') as error_file:
                output, error = result.communicate()
                error_file.write(str(error.strip().decode('utf-8')))

        else:
            # Adds the day to the list of checked dates
            checked_dates[unit].append(day)

            # Dumps the updated list of checked dates to a json file for use the next day
            with open(f'{auto_search_path}/checked_dates.json', 'w') as date_file:
                json.dump(checked_dates, date_file)

    return checked_dates


def main():
    results_path = '/home/jacob/search'
    auto_search_path = results_path + '/Autosearch'
    tl.make_path(auto_search_path)

    # Checks to see if the program is already running (maybe the dataset it's checking is quite large or long)
    try:
        with open(f'{auto_search_path}/pid.txt', 'r') as existing_pid_file:
            pid = int(existing_pid_file.readline())
            if psutil.pid_exists(pid):
                exit()
            else:
                raise FileNotFoundError

    # Runs the program normally if it isn't running already
    except FileNotFoundError:
        with open(f'{auto_search_path}/pid.txt', 'w') as pid_file:
            pid_file.write(str(os.getpid()))

        try:
            with open(f'{auto_search_path}/checked_dates.json', 'r') as date_file:
                checked_dates = json.load(date_file)

        # If the list ever gets deleted by accident or something
        except FileNotFoundError:
            checked_dates = {}

        # Running the main program on each of the detectors

        # THOR1
        checked_dates = search('THOR1', checked_dates, '/media/AllDetectorData/Detectors/THOR' + '/THOR1/Data',
                               results_path, auto_search_path)

        # THOR2
        checked_dates = search('THOR2', checked_dates, '/media/AllDetectorData/Detectors/THOR' + '/THOR2/Data',
                               results_path, auto_search_path)

        # THOR3
        checked_dates = search('THOR3', checked_dates, '/media/AllDetectorData/Detectors/THOR' + '/THOR3/Data',
                               results_path, auto_search_path)

        # THOR4
        checked_dates = search('THOR4', checked_dates, '/media/AllDetectorData/Detectors/THOR' + '/THOR4/Data',
                               results_path, auto_search_path)

        # THOR5
        checked_dates = search('THOR5', checked_dates, '/media/AllDetectorData/Detectors/THOR' + '/THOR5/Data',
                               results_path, auto_search_path)

        # THOR6
        checked_dates = search('THOR6', checked_dates, '/media/AllDetectorData/Detectors/THOR' + '/THOR6/Data',
                               results_path, auto_search_path)

        # GODOT
        checked_dates = search('GODOT', checked_dates, '/media/AllDetectorData/Detectors/SANTIS/Data',
                               results_path, auto_search_path)

        # SANTIS
        checked_dates = search('SANTIS', checked_dates, '/media/AllDetectorData/Detectors/GODOT/Data',
                               results_path, auto_search_path)

        # Deletes the pid file
        os.remove(f'{auto_search_path}/pid.txt')


if __name__ == '__main__':
    main()
    