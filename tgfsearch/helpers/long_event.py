"""A helper class used to keep track of and process long events."""
import numpy as np


class LongEvent:
    def __init__(self, start, length, z_scores, day_bins):
        self.start = int(start)
        self.stop = int(start + length - 1) if length > 1 else int(start + length)  # Non-inclusive (for slicing)
        self.length = int(length)

        glow_scores = z_scores[self.start:self.stop]
        self.peak_index = np.argmax(glow_scores) + self.start
        self.highest_score = glow_scores[self.peak_index - self.start]

        self.start_sec = day_bins[self.start]
        self.stop_sec = day_bins[self.stop]

        self.lm_files = {}

    # String casting magic method
    def __str__(self):
        return f'LongEvent({self.start}, {self.length}, ...)'

    # Debugging string magic method
    def __repr__(self):
        return self.__str__()
