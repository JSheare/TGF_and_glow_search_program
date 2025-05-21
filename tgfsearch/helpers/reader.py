"""A class that wraps and serves as an easy interface for the data reader."""

import tgfsearch.data_reader_annotated as dr


class Reader:
    """A class that wraps and serves as an easy interface for the data reader.

    Attributes
    ----------
    passtime : dict
        A dictionary containing information needed to import the subsequent list mode files properly (if applicable).

    """
    def __init__(self):
        # Data file documentation
        # We should probably move this to the data reader module in the future, but here is fine for now
        # Data files are broken down into a number of frames (or buffers)
        # A frame is a chunk of 2047 events (recorded in a single-line json-formatted string) preceded by either one or
        #   four time tags (old firmware uses one, more modern firmware uses four)
        # First tag: the time that the computer asked for the buffer
        # Second tag: the time that the computer received the "full buffer" message
        # Third tag: the time that the computer asked for the frame to be sent
        # Fourth tag: the time when the buffer finished arriving
        # Passtime entries:
        # lastsod: The second of day as calculated for the last event in the previous frame
        # ppssod: The second of day as calculated for the last GPS pulse per second of the previous frame
        # lastunix: Unix time (epoch seconds) for the last event of the previous frame (directly equivalent to lastsod
        #   regardless of data)
        # ppsunix: Unix time (epoch seconds) for the last GPS pulse per second of the previous frame (directly
        #   equivalent to ppssod, regardless of data)
        # lastwc: Bridgeport wall clock for the last event in the previous frame (no rollover corrections)
        # ppswc: Bridgeport wall clock for the last GPS pulse per second of the previous frame (no rollover corrections)
        # hz: Sampling rate of the analog to digital converter that records events
        # started: flag for whether or not there is a previous frame. If 0, current passtime values will be ignored
        # See reset() method for the default values of these entries
        self.passtime = {}
        self.reset()  # Sets passtime to the default

    def __call__(self, file_name, reset_after=False):
        """Defines the behavior for calling a class instance like a function."""
        return self.read(file_name, reset_after)

    def reset(self):
        """Resets the reader instance back to its default state."""
        self.passtime['lastsod'] = -1.0
        self.passtime['ppssod'] = -1.0
        self.passtime['lastunix'] = -1.0
        self.passtime['ppsunix'] = -1.0
        self.passtime['lastwc'] = -1
        self.passtime['ppswc'] = -1
        self.passtime['hz'] = 8e7  # Always this value for now
        self.passtime['started'] = 0

    def read(self, file_name, reset_after=False, clean_energy=False):
        """Reads the data file with the given name and returns the data as a pandas dataframe.

        Parameters
        ----------
        file_name : str
            The name of the data file to be read.

        reset_after : bool
            Optional. If True, resets the reader instance back to its default state after the data file has been
            read.
        clean_energy : bool
            Optional. If True, asks the data reader to strip out maximum and low energy counts.

        Returns
        -------
        pandas.core.frame.DataFrame
            A pandas dataframe containing the given file's data.

        """

        results = dr.fileNameToData(file_name, self.passtime, killcr=clean_energy)
        if type(results) is tuple:
            # List mode data
            data = results[0]
            self.passtime = results[1]
        else:
            # Trace data
            data = results

        if reset_after:
            self.reset()

        return data
