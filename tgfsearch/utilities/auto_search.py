"""A script that automatically runs searches on unexamined data."""
import os as os
import sys as sys
import json as json
import glob as glob
import psutil as psutil
import traceback as traceback

# Adds grandparent directory to sys.path. Necessary to make the imports below work when running this file as a script
grandparent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
if grandparent_dir not in sys.path:
    sys.path.append(grandparent_dir)

import tgfsearch.tools as tl
from tgfsearch.search import program


def search(unit, checked_dates, data_path, results_path, auto_search_path):
    if unit not in checked_dates:
        checked_dates[unit] = []

    queue = []
    days = glob.glob(f'{data_path}/*')
    days.sort()
    for day in days:
        if day[-6:].isnumeric():
            # If folder contains data and hasn't already been checked
            if len(os.listdir(day)) > 0 and day[-6:] not in checked_dates[unit]:
                queue.append(day[-6:])

    # Check the most recent stuff first
    queue = queue[::-1]

    for day in queue:
        mode_info = []
        if unit == 'THOR1':
            # Control flow for THOR1 NOAA flights
            if 220930 <= int(day) <= 230208:
                mode_info = ['--aircraft', '--skcali']

        elif unit == 'SANTIS':
            # Control flow for SANTIS instrument becoming the CROATIA instrument
            if int(day) > 211202:
                unit = 'CROATIA'

        # Runs search for the day
        print(f'Running search: {day} {unit}...')
        old_pwd = os.getcwd()
        os.chdir(results_path)
        try:
            program(day, day, unit, mode_info)
            os.chdir(old_pwd)

            # Adds the day to the list of checked dates
            checked_dates[unit].append(day)

            # Dumps the updated list of checked dates to a json file for use the next day
            with open(f'{auto_search_path}/checked_dates.json', 'w') as date_file:
                json.dump(checked_dates, date_file)

        except Exception as ex:
            # Writes any errors that the search runs into to text files so that they can be examined/fixed later
            print(f'Search encountered the following error: {ex}')
            tl.make_path(f'{auto_search_path}/Error Logs')
            with open(f'{auto_search_path}/Error Logs/{unit}_{day}_Error.txt', 'w') as error_file:
                error_file.write(traceback.format_exc())

    return checked_dates


def main():
    if len(sys.argv) >= 2:
        results_path = sys.argv[1].replace('\\', '/')
    else:
        results_path = os.getcwd().replace('\\', '/')

    if len(sys.argv) >= 3:
        auto_search_path = sys.argv[2].replace('\\', '/')
    else:
        auto_search_path = (os.getcwd() + '\\Autosearch').replace('\\', '/')

    if len(sys.argv) >= 4:
        data_dir = sys.argv[3].replace('\\', '/')
    else:
        data_dir = '/media/tgfdata/Detectors'

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

        # Running the program on each of the detectors

        # THOR1
        checked_dates = search('THOR1', checked_dates, data_dir + '/THOR/THOR1/Data',
                               results_path, auto_search_path)

        # THOR2
        checked_dates = search('THOR2', checked_dates, data_dir + '/THOR/THOR2/Data',
                               results_path, auto_search_path)

        # THOR3
        checked_dates = search('THOR3', checked_dates, data_dir + '/THOR/THOR3/Data',
                               results_path, auto_search_path)

        # THOR4
        checked_dates = search('THOR4', checked_dates, data_dir + '/THOR/THOR4/Data',
                               results_path, auto_search_path)

        # THOR5
        checked_dates = search('THOR5', checked_dates, data_dir + '/THOR/THOR5/Data',
                               results_path, auto_search_path)

        # THOR6
        checked_dates = search('THOR6', checked_dates, data_dir + '/THOR/THOR6/Data',
                               results_path, auto_search_path)

        # GODOT
        checked_dates = search('GODOT', checked_dates, data_dir + '/GODOT/Data',
                               results_path, auto_search_path)

        # SANTIS
        checked_dates = search('SANTIS', checked_dates, data_dir + '/SANTIS/Data',
                               results_path, auto_search_path)

        # Deletes the pid file
        os.remove(f'{auto_search_path}/pid.txt')


if __name__ == '__main__':
    main()
    