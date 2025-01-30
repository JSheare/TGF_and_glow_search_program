"""A class used to keep track of and process long events."""
import numpy as np


class LongEvent:
    """A class used to store all relevant information about a potential long event.

    Parameters
    ----------
    start : int
        The index of the histogram bin which corresponds to the start of the event.
    length : int
        The number of histogram bins that make up the event.
    z_scores : numpy.ndarray
        An array containing the z-scores for each bin in the daily histogram.
    day_bins : numpy.ndarray
        An array containing the bins for the daily histogram.

    Attributes
    ----------
    start : int
        The index of the histogram bin which corresponds to the start of the event.
    length : int
        The number of histogram bins that make up the event.
    stop : int
        The index of the histogram bin which corresponds to the end of the event.
    peak_index : int
        The location of the bin with the largest z-score among all the bins comprising the event.
    highest_score : float
        The largest z-score in the event.
    start_sec : int
        The beginning of the event in seconds since the beginning of the day.
    stop_sec : int
        The end of the event in seconds since the beginning of the day.
    lm_files : dict
        The list mode file(s) that the event is associated with.

    """

    def __init__(self, start, length, z_scores, day_bins):
        self.start = int(start)
        self.length = int(length)
        self.stop = int(start + length - 1) if self.length > 1 else int(start + length)
        self.peak_index, self.highest_score = self._highest_zscore(z_scores)
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
        peak_index = np.argmax(glow_scores) + self.start
        return peak_index, highest_score

    def _beginning_and_end_seconds(self, day_bins):
        """Calculates the beginning and end of an event in seconds."""
        glow_times = day_bins[self.start:self.stop]
        first_sec = glow_times[0]
        length = self.length * int(day_bins[1] - day_bins[0])
        last_sec = first_sec + length
        return first_sec, last_sec
