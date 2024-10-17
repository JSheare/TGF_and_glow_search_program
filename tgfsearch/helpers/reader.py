"""A class that wraps and serves as an easy interface for the data reader."""

import pandas as pd

import tgfsearch.DataReaderTimetrack2 as Dr


class Reader:
    """A class that wraps and serves as an easy interface for the data reader.

    Attributes
    ----------
    passtime : dict
        A dictionary containing information needed to import the subsequent list mode files properly (if applicable).

    """
    def __init__(self):
        self.passtime = {}
        self.reset()  # Sets passtime to the default

    def __call__(self, file_name, reset_after=False):
        """Defines the behavior for calling a class instance like a function."""
        return self.read(file_name, reset_after)

    def reset(self):
        """Resets the reader instance back to its default state."""
        self.passtime = {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0, 'ppsunix': -1.0, 'lastwc': -1, 'ppswc': -1,
                         'hz': 8e7, 'started': 0}

    def read(self, file_name, reset_after=False):
        """Reads the data file with the given name and returns the data as a pandas dataframe.

        Parameters
        ----------
        file_name : str
            The name of the data file to be read.

        reset_after : bool
            Optional. If True, resets the reader instance back to its default state after the data file has been
            read.

        Returns
        -------
        pandas.core.frame.DataFrame
            A pandas dataframe containing the given file's data.

        """

        data, self.passtime = Dr.fileNameToData(file_name, self.passtime)
        if reset_after:
            self.reset()

        return data
