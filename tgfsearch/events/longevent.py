"""A class used to keep track of and process long events."""
import numpy as np

import tgfsearch.parameters as params


class LongEvent:
    """Object used to store all relevant information about long events.

    Parameters
    ----------
    start : int
        The index of the histogram bin which corresponds to the start of the event.
    length : int
        The number of histogram bins which make up an event.
    z_scores : list
        A list containing the z-scores for each bin in the daily histogram.
    day_bins : np.array
        An array containing the bins for the daily histogram.

    Attributes
    ----------
    start : int
        The index of the histogram bin which corresponds to the start of the event.
    length : int
        The number of histogram bins which make up an event.
    stop : int
        The index of the histogram bin which corresponds to the end of the event.
    peak_index : int
        The location of the bin with the largest z-score among all the bins comprising the event.
    highest_score : float
        The largest z-score in the event.
    start_sec : int
        The beginning of the event in seconds.
    stop_sec : int
        The end of the event in seconds.
    lm_files : dict
        The list mode file(s) that the event is associated with.

    """

    def __init__(self, start, length, z_scores, day_bins):
        self.start = int(start)
        self.length = int(length)
        self.stop = int(start + length - 1) if self.length > 1 else int(start + length)
        self.peak_index = 0
        self.highest_score = self._highest_zscore(z_scores)
        self.start_sec, self.stop_sec = self._beginning_and_end_seconds(day_bins)
        self.lm_files = {}

    # String casting magic method
    def __str__(self):
        return f'Long event at index {self.start}'

    # Debugging string magic method
    def __repr__(self):
        return f'Long event; start:{self.start}; stop:{self.stop}; length:{self.length}'

    def _highest_zscore(self, z_scores):
        """Identifies the highest z-score and its corresponding bin for the event."""
        glow_scores = z_scores[self.start:self.stop]
        highest_score = np.max(glow_scores)
        self.peak_index = np.argmax(glow_scores) + self.start
        return highest_score

    def _beginning_and_end_seconds(self, day_bins):
        """Calculates the beginning and end of an event in seconds."""
        glow_times = day_bins[self.start:self.stop]
        first_sec = glow_times[0]
        length = self.length * params.BIN_SIZE
        last_sec = first_sec + length
        return first_sec, last_sec

    def find_lm_filenames(self, filelist_dict, extrema_dict):
        """Gets the names of the files that the event occurred in. Note: this code assumes that events are in
        chronological order.

        Parameters
        ----------
        filelist_dict : dict
            A dictionary containing a list of files for each scintillator that contributed to the event.
        extrema_dict : dict
            A dictionary containing a list of file time extrema for each scintillator that contributed to the event.

        """

        for scintillator in filelist_dict:
            if scintillator not in self.lm_files:
                event_file = ''
                files_examined = 0
                filetime_extrema = extrema_dict[scintillator]
                filelist = filelist_dict[scintillator]
                for i in range(len(filetime_extrema)):
                    first_time = filetime_extrema[i][0]
                    last_time = filetime_extrema[i][1]
                    if (first_time <= self.start_sec <= last_time) or (
                            first_time <= self.stop_sec <= last_time):
                        event_file = filelist[i]
                        # Warning: modifies the actual dictionaries because python passes dictionaries by reference
                        filelist_dict[scintillator] = filelist[files_examined:]
                        extrema_dict[scintillator] = filetime_extrema[files_examined:]
                        break

                    files_examined += 1

                self.lm_files[scintillator] = event_file
