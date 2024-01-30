import os as os
import subprocess as subprocess
import json as json
import glob as glob
import psutil as psutil
import sys as sys
search_directory = '/home/jacob/search/'  # This would need to be updated if search.py or its modules are moved
autosearch_directory = os.path.dirname(os.path.abspath(sys.argv[0])) + '/'
sys.path.insert(0, search_directory)
import search_module as sm


def search(detector, dates, path):
    queue = []
    days = glob.glob(f'{path}/*')
    days.sort()
    already_checked = dates[detector]
    for day in days:
        if day[-6:].isnumeric():
            has_data = True if len(os.listdir(day)) > 0 else False
            if has_data and day[-6:] not in already_checked:
                queue.append(day[-6:])

    # Check the most recent stuff first
    queue = queue[::-1]

    for day in queue:
        # Note: these strings will have to be updated if search.py is moved to another directory
        if detector == 'THOR1':
            # Control flow for THOR1 NOAA flights
            if 220930 <= int(day) <= 230208:
                command = ['python3', 'search.py', day, day, 'THOR1', 'aircraft', 'skcali']
            else:
                command = ['python3', 'search.py', day, day, 'THOR1']

        elif detector == 'SANTIS':
            # Control flow for SANTIS instrument becoming the CROATIA instrument
            if int(day) > 211202:
                command = ['python3', 'search.py', day, day, 'CROATIA']
            else:
                command = ['python3', 'search.py', day, day, 'SANTIS']

        else:
            command = ['python3', 'search.py', day, day, f'{detector}']

        # Runs search.py for the day
        old_pwd = os.getcwd()
        os.chdir(search_directory)
        result = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Uncomment this and comment the above line to get search.py console output (for console debugging)
        # result = subprocess.Popen(command)
        os.chdir(old_pwd)

        if result.wait() != 0:
            # Writes any errors that search.py runs into to text files so that they can be examined/fixed later
            if not os.path.exists(autosearch_directory + 'Error Logs'):
                os.makedirs(autosearch_directory + 'Error Logs')

            with open(autosearch_directory + f'Error Logs/{detector}_{day}_Error.txt', 'w') as error_file:
                output, error = result.communicate()
                error_file.write(str(error.strip().decode('utf-8')))

        else:
            # Adds the day to the list of checked dates
            dates[detector].append(day)

            # Dumps the updated list of checked dates to a json file for use the next day
            with open(autosearch_directory + 'checked_dates.json', 'w') as date_file:
                json.dump(dates, date_file)

    return dates


def main():
    # Checks to see if the program is already running (maybe the dataset it's checking is quite large or long)
    try:
        with open(autosearch_directory + 'pid.txt', 'r') as existing_pid_file:
            pid = int(existing_pid_file.readline())
            if psutil.pid_exists(pid):
                exit()
            else:
                raise FileNotFoundError

    # Runs the program normally if it isn't running already
    except FileNotFoundError:
        with open(autosearch_directory + 'pid.txt', 'w') as pid_file:
            pid_file.write(str(os.getpid()))

        try:
            date_file = open(autosearch_directory + 'checked_dates.json')
            checked_dates = json.load(date_file)
            date_file.close()
        # If the list ever gets deleted by accident or something
        # (mostly here because I don't want to put the list together by hand)
        except FileNotFoundError:
            checked_dates = {'THOR1': [], 'THOR2': [], 'THOR3': [],
                             'THOR4': [], 'THOR5': [], 'THOR6': [], 'GODOT': [], 'SANTIS': []}

        # Running the main program on each of the detectors

        # THOR1
        checked_dates = search('THOR1', checked_dates, sm.T_raw_data_loc() + '/THOR1/Data')

        # THOR2
        checked_dates = search('THOR2', checked_dates, sm.T_raw_data_loc() + '/THOR2/Data')

        # THOR3
        checked_dates = search('THOR3', checked_dates, sm.T_raw_data_loc() + '/THOR3/Data')

        # THOR4
        checked_dates = search('THOR4', checked_dates, sm.T_raw_data_loc() + '/THOR4/Data')

        # THOR5
        checked_dates = search('THOR5', checked_dates, sm.T_raw_data_loc() + '/THOR5/Data')

        # THOR6
        checked_dates = search('THOR6', checked_dates, sm.T_raw_data_loc() + '/THOR6/Data')

        # GODOT
        checked_dates = search('GODOT', checked_dates, sm.G_raw_data_loc())

        # SANTIS
        checked_dates = search('SANTIS', checked_dates, sm.S_raw_data_loc())

        # Deletes the pid file
        os.remove(autosearch_directory + 'pid.txt')


if __name__ == '__main__':
    main()
    