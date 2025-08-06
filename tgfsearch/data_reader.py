"""A module containing functions for reading UCSC TGF group data files.

Changelog:
Version 2/26/23
    - Changed sign of unixtimecorrection to subtraction in processatatiming()
    - Changed location of best time tag string in thor data lines[i-4] (first tag) in thorFileToData()
    - Added ucsc_timestring_to_datetime() routine to correctly pad microseconds in UCSC time tags with 0s on the left,
      needed by datetime.strptime()
    - Added multifilesToData() routine to read multiple files in a single dataset
    - Added warnings about unexpected behavior across frame boundaries to thorFileToData()
Version 5/6/24
    - Added handling of the case where there are columns for specific flags (like 'pps') instead of a 'flags' column
      like in normal THOR data in getdatafromlmthor()
    - Changed "if mode is" to "if mode ==" on the advice of new warnings.
Version 7/5/24
    - Fixed condition where there is a rollover between the last PPS of the prior frame and the end of the prior frame,
      which was messing up all the counts up to the first pps of the new frame.
Version 3/3/25
    - Added "killcr" feature for noisy detectors in France, kills overflow and unrealistically low-energy counts
    - Also made lmFileToData take and pass passtime consistently (don't know if needed)
Version 6/12/25
    - Changed architecture to use parsers for every data file format
    - Made SecondsOfDay consistent across day boundaries (resets to zero)
Version 7/9/25
    - Added get_passtime function, which returns a fresh passtime dictionary
    - Replaced all print statements with warnings, which can be ignored, caught, or redirected

"""
import gzip as gzip
import json as json
import numpy as np
import pandas as pd
import re as re
import warnings as warnings
from datetime import datetime


def get_passtime():
    """Returns a fresh instance of the passtime dictionary needed for reading consecutive frames/files.

    Passtime entries:
    - lastsod: The second of day as calculated for the last event in the previous frame
    - ppssod: The second of day as calculated for the last GPS pulse per second of the previous frame
    - lastunix: Unix time (epoch seconds) for the last event of the previous frame (directly equivalent to lastsod
        regardless of data)
    - ppsunix: Unix time (epoch seconds) for the last GPS pulse per second of the previous frame (directly equivalent to
      ppssod, regardless of data)
    - lastwc: Bridgeport wall clock for the last event in the previous frame (no rollover corrections)
    - ppswc: Bridgeport wall clock for the last GPS pulse per second of the previous frame (no rollover corrections)
    - frlen: the length of the first frame read
    - prevfr: The most recently parsed frame. Used to check for duplicate frame data
    - started: flag for whether there is a previous frame. If 0, current passtime values will be ignored

    Returns
    -------
    dict
        The passtime dictionary.

    """

    return {'lastsod': -1, 'ppssod': -1, 'lastunix': -1, 'ppsunix': -1, 'lastwc': -1, 'ppswc': -1, 'frlen': -1,
            'prevfr': None, 'started': 0}


def read_files(file_file):
    with open(file_file, 'r') as file:
        lines = file.read().splitlines()

    file_frames = []
    passtime = get_passtime()
    for file in lines:
        print(f'Reading file: {file}')
        data, passtime = read_file(file, passtime)
        file_frames.append(data)

    return pd.concat(file_frames, axis=0)


def read_file(file_name, passtime=None, killcr=0):
    if file_name[-2:] == 'gz':
        # Decompresses the file if it's compressed
        with gzip.open(file_name, 'rt') as file:
            lines = file.read().splitlines()
    else:
        with open(file_name, 'r') as file:
            lines = file.read().splitlines()

    if passtime is None:
        passtime = get_passtime()

    # This means it's a trace file, so go to the trace reading function
    if 'xtr' in file_name:
        return read_trace_file(lines), passtime

    # Only the json-formatted list mode files start with a 3 character (including a space) line
    if len(lines[0]) == 3:
        # Only the NRL software has a named pps column
        if 'pps' in lines[5]:
            parser = nrl_lm_parser(passtime, lines)
        else:
            parser = thor_lm_parser(passtime, lines)
    else:
        # Early list mode files alternate (Time) - (Buffer) - (Time) - (Buffer); the newer version has 4 time tags
        # before a buffer. So, we can determine the type by testing whether the first line is a time tag
        if len(lines[2]) < 50:
            parser = generic_lm_parser_new(passtime, lines)
        else:
            parser = generic_lm_parser_old(passtime, lines)

    return read_lm_file(passtime, parser, killcr), passtime


def ucsc_timestring_to_datetime(time_string):
    parts = time_string.split()
    parts[-1] = parts[-1].rjust(6, '0')
    return datetime.strptime(' '.join(parts), '%Y %m %d %H %M %S %f')


def custom_warning_format(message, category, filename, lineno, line=None):
    return f'Warning: {message}\n'


# Raises a warning with a custom message
def raise_warning(message):
    original = warnings.formatwarning
    try:
        # Setting a custom format for the warning string
        warnings.formatwarning = custom_warning_format
        warnings.warn(message)
    finally:
        warnings.formatwarning = original


# A generator that parses old-style generic list mode data files (previously mode 1)
def generic_lm_parser_old(passtime, lines):
    """Format Info:
    - Every file contains a single blank line, a file creation timestamp line (?), and then all the frame blocks.
    - Each frame block consists of a frame data line and a timestamp line.
    - The timestamp line contains a timestamp corresponding to the time when the frame was requested by the computer (?).
    - The frame data line is just a collection of space-separated words representing the frame's data.
    - The first two words are the detector's serial number and the number of events in the frame.
    - The remaining words are the actual data.
    - Each event consists of three words. The first is the event energy, and the last two are 16-bit halves of the
        32-bit integer representing the event's wall clock tick.
    - The wall clock is either reset after every frame or rolls over after 2^32 ticks, depending on the configuration.
    - The ADC (theoretically) samples at a frequency of 8*10^7 hz, so this is the rate that the wall clock ticks at.
    """

    # Starting on line three (index two). Each frame is two lines long
    for i in range(2, len(lines), 2):
        frame_header = ucsc_timestring_to_datetime(lines[i + 1])  # Header line is one ahead of the lm string line
        # Splitting the lm string into individual values. Throwing out the first two, which are unneeded
        frame_values = [int(x) for x in lines[i].split(' ')[2:]]
        # Making each row. Every row is three values long
        frame_data = pd.DataFrame(np.reshape(frame_values, (int(len(frame_values)/3), 3)))
        # Recreating full 32-bit wall clock values from two 16-bit words (one fine, one coarse)
        # 65536 is 2^16. Multiplying an int by it is the same as a 16 bit rightward shift
        wc = frame_data[1] + frame_data[2] * 65536
        energy = frame_data[0]
        flags = frame_data[0] * 0
        frame_data = pd.concat([energy.rename('energy'), wc.rename('wc'), flags.rename('flags')], axis=1)
        # Eliminating duplicate data from the previous frame, if it's present
        if passtime['prevfr'] is not None and frame_data['wc'][0] == passtime['prevfr']['wc'][0]:
            frame_data.drop([j for j in range(0, len(passtime['prevfr'].index))], inplace=True)
            frame_data.reset_index(inplace=True)

        passtime['prevfr'] = frame_data

        # Yielding: the frame's time header, the frame's data, the adc sampling rate, and the rollover period
        yield frame_header, frame_data, 8e7, 2**32


# A generator that parses new-style generic list mode data files (previously mode 2)
def generic_lm_parser_new(passtime, lines):
    """Format Info:
    - Every file contains two lines of metadata followed by all the frame blocks.
    - Each frame block consists of six lines: a buffer line, four timestamp lines, and a frame data line.
    - The buffer line is just two numbers (separated by a space) that record which data buffer the frame was read from
        (first number) and whether that data buffer was full (second number).
    - The four timestamp lines record (in order):
        - When the frame was requested by the computer.
        - When the computer received the "full buffer" message.
        - When the computer asked for the frame to be sent.
        - When the frame finished arriving.
    - The frame data line is just a collection of space-separated words representing the frame's data.
    - The first three words are (something), the detector's serial number, and the number of events in the frame.
    - The remaining words are actual data.
    - Note: there is always one dummy event at the beginning of the frame that should be ignored.
    - Each event is six words long:
        - (something).
        - The event's energy.
        - One third (16-bits) of the 48-bit integer representing the wall clock tick (fine).
        - Another third (16-bits) of the 48-bit integer representing the wall clock tick (coarse).
        - The last third (16-bits) of the 48-bit integer representing the wall clock tick (very coarse).
        - The event's flags.
    - The wall clock rolls over after 2^48 ticks.
    - The ADC (theoretically) samples at a frequency of 8*10^7 hz, so this is the rate that the wall clock ticks at.
    """

    rollover_period = 2**48
    # Starting on line eight (index seven). Each frame is six lines long
    for i in range(7, len(lines), 6):
        frame_header = ucsc_timestring_to_datetime(lines[i - 4])  # Header line is four behind the lm string line
        # Splitting the lm string into individual values. Throwing out the first nine, which are unneeded
        frame_values = [int(x) for x in lines[i].split(' ')[9:]]
        # Making each row. Every row is six values long
        frame_data = pd.DataFrame(np.reshape(frame_values, (int(len(frame_values) / 6), 6)))
        # Recreating full 64-bit wall clock values from three 16-bit words (one fine, one coarse, and one very coarse)
        # Casting to a 64-bit integer is necessary to prevent integer overflow
        # 65536 is 2^16. Multiplying an int by it is the same as a 16 bit rightward shift
        wc = (frame_data[2].astype('int64') +
              frame_data[3].astype('int64') * 65536 +
              frame_data[4].astype('int64') * 65536**2)  # Shifting right by 32 bits
        energy = frame_data[1]
        # Format string converts x to its string representation as a binary integer with at least eight bits
        pps = pd.Series([int('{0:08b}'.format(x)[-8]) for x in frame_data[5]])
        flags = frame_data[5]
        frame_data = pd.concat([energy.rename('energy'), wc.rename('wc'), pps.rename('PPS'), flags.rename('flags')],
                               axis=1)
        frame_data['gpsSync'] = False
        # Eliminating duplicate data from the previous frame. This usually happens when buffers that are only partially
        # full are read out
        if passtime['prevfr'] is not None:
            # Data is duplicated from some midpoint to the end of the frame
            diff = np.where(frame_data['wc'].diff() < 0)[0]
            # Checking that this isn't a rollover, and dropping the duplicate rows if it isn't
            if len(diff) > 0 and frame_data['wc'][diff[0]] - frame_data['wc'][diff[0] - 1] < rollover_period/4:
                frame_data.drop([j for j in range(diff[0], len(frame_data.index))], inplace=True)
                frame_data.reset_index(inplace=True)

            # The whole frame is duplicated data from before the previous frame
            if passtime['prevfr']['wc'][0] > frame_data['wc'][len(frame_data.index) - 1]:
                # Checking that this isn't a rollover, and skipping the whole frame if it isn't
                if passtime['prevfr']['wc'][0] - frame_data['wc'][len(frame_data.index) - 1] < rollover_period/4:
                    continue

        passtime['prevfr'] = frame_data

        # Yielding: the frame's time header, the frame's data, the adc sampling rate, and the rollover period
        yield frame_header, frame_data, 8e7, rollover_period


# A generator that parses Thor list mode data files (previously mode 0)
def thor_lm_parser(passtime, lines):
    """Format Info:
    - Every file is just a collection of frame blocks.
    - Each frame block consists of six lines: a buffer line, four timestamp lines, and a frame data line.
    - The buffer line is just two numbers (separated by a space) that record which data buffer the frame was read from
        (first number) and whether that data buffer was full (second number).
    - The four timestamp lines record (in order):
        - When the frame was requested by the computer.
        - When the computer received the "full buffer" message.
        - When the computer asked for the frame to be sent.
        - When the frame finished arriving.
    - The frame data line is a string containing all the frame's data.
    - After the detector's serial number and a space, the rest of the string is just the json-formatted data.
    - The wall clock rolls over after 2^36 ticks.
    - The ADC (theoretically) samples at a frequency of 8*10^7 hz, so this is the rate that the wall clock ticks at.
    """

    # Starting on line six (index five). Each frame is six lines long
    for i in range(5, len(lines), 6):
        frame_header = ucsc_timestring_to_datetime(lines[i - 4])  # Header line is four behind the lm string line
        frame_data = pd.DataFrame(json.loads(re.sub('eRC[0-9]{4} ', '', lines[i]))['lm_data'])
        # Format string converts x to its string representation as a binary integer with at least eight bits
        frame_data['PPS'] = [int('{0:08b}'.format(x)[-8]) for x in frame_data['flags']]
        frame_data['gpsSync'] = False
        # Yielding: the frame's time header, the frame's data, the adc sampling rate, and the rollover period
        yield frame_header, frame_data, 8e7, 2**36


# A generator that parses NRL json-style list mode data files (previously mode 3)
def nrl_lm_parser(passtime, lines):
    """Format Info:
    - Every file is just a collection of frame blocks.
    - Each frame block consists of six lines: a buffer line, four timestamp lines, and a frame data line.
    - The buffer line is just two numbers (separated by a space) that record which data buffer the frame was read from
        (first number) and whether that data buffer was full (second number).
    - The four timestamp lines record (in order):
        - When the frame was requested by the computer.
        - When the computer received the "full buffer" message.
        - When the computer asked for the frame to be sent.
        - When the frame finished arriving.
    - The frame data line is a string containing all the frame's data.
    - After the detector's serial number and a space, the rest of the string is just the json-formatted data.
    - Note: unlike Thor data, this format has a dedicated PPS (pulse per second) column built in.
    - The wall clock rolls over after 2^48 ticks.
    - The ADC (theoretically) samples at a frequency of 8*10^7 hz, so this is the rate that the wall clock ticks at.
    """

    # Starting on line six (index five). Each frame is six lines long
    for i in range(5, len(lines), 6):
        frame_header = ucsc_timestring_to_datetime(lines[i - 4])  # Header line is four behind the lm string line
        frame_data = pd.DataFrame(json.loads(re.sub('eRC[0-9]{4} ', '', lines[i]))['lm_data'])
        frame_data.rename(columns={'pps': 'PPS'}, inplace=True)
        frame_data['gpsSync'] = False
        # Yielding: the frame's time header, the frame's data, the adc sampling rate, and the rollover period
        yield frame_header, frame_data, 8e7, 2**48


# Reads a list mode file and returns its data
def read_lm_file(passtime, parser, killcr):
    data_list = []
    started = passtime['started']
    for header, data, adc_rate, rollover_period in parser:
        last_sod = passtime['lastsod']
        process_data_timing(passtime, header, data, adc_rate, rollover_period)
        first_count = data['SecondsOfDay'][0]
        if started:
            dt = first_count - last_sod
            if dt > 0.5:
                raise_warning(f'long gap ({dt}s) between frames at {header}')

            if dt < 0:
                raise_warning(f'anomalous clock backwards at {header}')

        if 'energy' in data.columns:
            data.rename(columns={'energy': 'energies'}, inplace=True)

        if killcr:
            data = data[data['energies'] < 65000]
            data = data[data['energies'] > 100]

        data_list.append(data)
        started = passtime['started']

    return pd.concat(data_list)


# Processes the timing information for a single frame of list mode data
def process_data_timing(passtime, header, data, adc_rate, rollover_period):
    last_index = len(data.index) - 1

    # Every time we encounter the wallclock going significantly backwards (1/4 or range), assume it's a rollover.
    # Count how many times it's rolled over and multiply that by the rollover period to get each row's correction
    # This is about making the wc column monotonically increasing
    rollover_correction = (data['wc'].diff() < -rollover_period/4).cumsum() * rollover_period
    data['wc'] += rollover_correction

    # Get the header timestamp in unix time, and unix time at start of second, and start of day
    header_unix = header.timestamp()
    # Subtracts out 00:00:00 of the earliest day in UNIX time.
    daystart_unix = datetime(header.year, header.month, header.day).timestamp()

    # Compensating for a missing frame right before the current one
    # If there's a gap of more than 1/4 the length of an average frame between the current frame and the last frame,
    # assume that the previous frame is missing
    if passtime['started']:
        first_count = header_unix - (data['wc'][last_index] - data['wc'][0])/adc_rate
        if first_count - passtime['lastunix'] > passtime['frlen']/4:
            raise_warning(f'previous frame may be missing at {header}')
            passtime['started'] = 0

    # Checking for a wallclock rollover before the start of the frame
    if passtime['started'] and data['wc'][0] - passtime['ppswc'] < -rollover_period/4:
        # Checking for a rollover between the last frame's pps and the final count
        if passtime['lastwc'] - passtime['ppswc'] < -rollover_period/4:
            passtime['lastwc'] += rollover_period

        data['wc'] += rollover_period

    # Checking for a wallclock reset between frames (older list mode data can do this)
    if passtime['started'] and data['wc'][0] < passtime['lastwc']:
        passtime['started'] = 0

    # Get PPS events (if available)
    if 'PPS' in data:
        pps = np.where(data['PPS'])[0]
        n_pps = len(pps)
    else:
        pps = None
        n_pps = 0

    # First method: in-frame only, last count = time tag. Interpolate all data in the frame according to the assumption
    # that the last event actually happened at header_unix seconds
    m1_unix = header_unix - (data['wc'][last_index] - data['wc']) / adc_rate

    # Second method: across-frame, not using PPS. Interpolate all data in the frame using the time of the last event in
    # the previous frame as a starting point.
    if passtime['started']:
        m2_unix = passtime['lastunix'] + (data['wc'] - passtime['lastwc']) / adc_rate
    else:
        m2_unix = m1_unix * 0 - 1  # Column full of -1s

    # If it's the very first buffer (or there aren't enough pps points to do an intra-frame interpolation with), base
    # everything on the header. The last event in the frame is assumed to have happened at exactly the header time.
    # The wallclock rate is assumed to be precise
    if n_pps < 2:
        if not passtime['started']:
            data['UnixTime'] = m1_unix
        else:
            data['UnixTime'] = m2_unix

        if 'PPS' in data:
            raise_warning(f'frame with {n_pps} PPS at {header}')

    # Otherwise, use pps to get a more accurate times
    else:
        # Finding the frequencies in between each pair of PPS

        # This column keeps track of them. Using the adc sampling rate as the default value
        frequencies = m1_unix * 0 + adc_rate
        drift_err = 1e4

        # Intra-frame pps interpolation
        for i in range(n_pps - 1):
            # The number of wallclock ticks between two consecutive pps points
            test_hz = data['wc'][pps[i + 1]] - data['wc'][pps[i]]

            # Checking to see that the actual wallclock tick frequency is within a certain acceptable error
            if abs(test_hz - adc_rate) < drift_err:
                # Updating the frequency record column for these two pps points
                frequencies[pps[i]:pps[i + 1] + 1] = test_hz
                data.loc[pps[i]:pps[i + 1] + 1, 'gpsSync'] = True
            else:
                # If not synced, try guessing that one PPS was skipped. Otherwise, out of sync.
                if abs(test_hz / 2 - adc_rate) < drift_err:
                    frequencies[pps[i]:pps[i + 1] + 1] = test_hz / 2
                    data.loc[pps[i]:pps[i + 1] + 1, 'gpsSync'] = True
                    raise_warning(f'missing one PPS at {header}. Freq: {test_hz}; Freq/2: {test_hz / 2}')
                else:
                    data.loc[pps[i]:pps[i + 1] + 1, 'gpsSync'] = False
                    raise_warning(f'GPS sync failed at {header}. Factor: {test_hz / adc_rate}')

            if i == 0:
                frequencies[0:pps[0]] = test_hz

        # Inter-frame pps interpolation
        # If prior PPS is available from before the frame, improve the starting events before first PPS:
        if passtime['started']:
            test_hz = data['wc'][pps[0]] - passtime['ppswc']
            if abs(test_hz - adc_rate) < drift_err:
                frequencies[0:pps[0]] = test_hz
                data.loc[0:pps[0], 'gpsSync'] = True
            else:
                # If not synced, try guessing that one PPS was skipped. Otherwise, out of sync.
                if abs(test_hz / 2 - adc_rate) < drift_err:
                    frequencies[0:pps[0]] = test_hz / 2
                    data.loc[0:pps[0], 'gpsSync'] = True
                    raise_warning(f'missing one PPS across frame boundary at {header}. '
                                  f'Freq: {test_hz}; Freq/2: {test_hz / 2}')
                else:
                    data.loc[0:pps[0], 'gpsSync'] = False
                    raise_warning(f'GPS sync failed at {header}. Factor: {test_hz / adc_rate}')

        # Extrapolating that events after the last pps have the same frequency as those right before
        frequencies[pps[-1]:] = frequencies[pps[-1] - 1]

        # Now get the full array of unix time by two methods.
        # Start with the reverse method (the only one when it's the very first buffer). We start at the last pps and
        # walk backwards, applying the interpolated wallclock frequency corrections
        index = 0  # index==0 is the *last* pps
        m3_unix = m1_unix * 0
        rev_pps = np.flip(pps)  # rev_pps[index] = i, as used below
        for i in rev_pps:
            if index == 0:  # Last pps in frame, first to be examined
                m3_unix[i] = round(m1_unix[i])
                m3_unix[i:] = m3_unix[i] + (data['wc'][i:] - data['wc'][i]) / frequencies[i:]
            # All but the very first pps
            else:
                known = rev_pps[index - 1]
                m3_unix[i:known] = m3_unix[known] - (data['wc'][known] - data['wc'][i:known]) / frequencies[i:known]

            # First pps in frame, last to be examined
            if index == n_pps - 1:
                m3_unix[0:i] = m3_unix[i] - (data['wc'][i] - data['wc'][0:i]) / frequencies[0:i]

            index += 1

        # Forward method if available. We start at the first pps (retrieved from passtime) and walk forwards, applying
        # the interpolated wallclock frequency corrections
        m4_unix = m1_unix * 0 - 1
        if passtime['started']:
            index = 0
            for i in pps:
                m4_unix[i] = round(m2_unix[i])
                # First pps in frame, first to be examined
                if index == 0:
                    m4_unix[0:i + 1] = (passtime['ppsunix'] +
                                        (data['wc'][0:i + 1] - passtime['ppswc']) / frequencies[0:i + 1])
                    # m4_unix[0:i] = m4_unix[i] - (data['wc'][i] - data['wc'][0:i])/frequencies[0:i]

                # All but the very last pps
                if index < len(pps) - 1:
                    m4_unix[i:pps[index + 1]] = (m4_unix[i] + (data['wc'][i:pps[index + 1]] - data['wc'][i]) /
                                                 frequencies[i:pps[index + 1]])
                # Last pps in frame, last to be examined
                else:
                    m4_unix[i:] = m4_unix[i] + (data['wc'][i:] - data['wc'][i]) / frequencies[i:]

                index = index + 1

        if passtime['started']:
            data['UnixTime'] = m4_unix
        else:
            data['UnixTime'] = m3_unix

    # Updating dataframe
    data['SecondsOfDay'] = data['UnixTime'] - daystart_unix
    # Checking to see if we've crossed a day boundary in this frame. When this happens, SecondsOfDay goes partially
    # negative. Counts in the current day become count - 86400
    negatives = np.where(data['SecondsOfDay'] < 0)[0]
    if len(negatives) > 0 and data['SecondsOfDay'].max() < 5000:
        # Updating data so that when the day boundary is crossed, SecondsOfDay will start over again from zero
        data.iloc[0:negatives[-1] + 1, data.columns.get_loc('SecondsOfDay')] += 86400

    # Undoing rollover correction
    data['wc'] -= rollover_correction

    # Updating passtime
    passtime['lastunix'] = data['UnixTime'][last_index]
    passtime['lastsod'] = data['SecondsOfDay'][last_index]
    passtime['lastwc'] = data['wc'][last_index]
    if n_pps > 0:
        passtime['ppswc'] = data['wc'][pps[-1]]
        passtime['ppsunix'] = data['UnixTime'][pps[-1]]
        passtime['ppssod'] = data['SecondsOfDay'][pps[-1]]
    else:
        passtime['ppswc'] = -1
        passtime['ppsunix'] = -1
        passtime['ppssod'] = -1

    # Recording the length of the very first frame read. Used to check for prior missing frames later
    if not passtime['started']:
        passtime['frlen'] = data['SecondsOfDay'][last_index] - data['SecondsOfDay'][0]

    passtime['started'] = 1


# Reads a trace file and returns its data
def read_trace_file(lines):
    data_list = []
    for i in range(4, len(lines), 5):
        data = pd.DataFrame.from_dict(json.loads(re.sub("eRC[0-9]{4} [0-9]", "", lines[i])))
        # There should be some significance to how this affects the time relationships?
        data['BufferNo'] = int(lines[i][8:9])
        data['DateTime'] = datetime.strptime(lines[i - 1], "%Y %m %d %H %M %S %f")
        data['Seconds'] = [x * 1.25e-8 for x in range(len(data))]
        data_list.append(data)

    return pd.concat(data_list)
