"""A class used to keep track of and process short events."""
import numpy as np

import tgfsearch.parameters as params
import tgfsearch.tools as tl


class ShortEvent:
    """A class used to store all relevant information about a potential short event.

    Parameters
    ----------
    event_start : int
        The index of the time array which corresponds to the start of the event.
    event_length : int
        The number of entries in the time array which make up the event.
    scintillator : str
        A string corresponding to the scintillator which the event was found in.

    Attributes
    ----------
    start : int
        The index of the time array which corresponds to the start of the event.
    length : int
        The number of entries in the time array which make up the event.
    stop : int
        The index of the time array which corresponds to the end of the event.
    scintillator : str
        A string corresponding to the scintillator which the event was found in.
    lm_files : dict
        The list mode file(s) associated with the event.
    traces
        The time id(s) of the trace(s) associated with the event.
    len_subscore : int
        The event's length score given by the ranking system.
    clumpiness_subscore : int
        The event's clumpiness score given by the ranking system.
    hel_subscore : int
        The event's high energy leading count score given by the ranking system.
    weather_subscore : int
        The event's weather score given by the ranking system.
    total_score : int
        The event's total score as determined by the subscores. Calculated by the ranking system.
    number : int
        The event's number in relation to the others (chronologically).
    rank : int
        The event's rank among all the other events found.

    """

    def __init__(self, event_start, event_length, scintillator):
        self.start = int(event_start)
        self.length = int(event_length)
        self.stop = int(event_start + event_length)
        self.scintillator = scintillator
        self.lm_files = {}
        self.traces = {}

        # Ranking scores
        self.len_subscore = 0
        self.clumpiness_subscore = 0
        self.hel_subscore = 0
        self.weather_subscore = 0
        self.total_score = 0

        self.number = 0
        self.rank = 0

    # String casting overload
    def __str__(self):
        return f'{self.scintillator} short event at index {self.start}'

    # Debugging string dunder
    def __repr__(self):
        return f'{self.scintillator} short event; start:{self.start}; stop:{self.stop}; length:{self.length}'

    def _calculate_subscores(self, detector, weather_cache, times, energies):
        """Calculates the event's ranking subscores and updates them."""
        clumpiness = 0
        clump_counts = 0

        high_energy_lead = 0
        leading_counts = 0

        for i in range(self.start, self.stop):
            if i > self.start:
                difference = times[i] - times[i - 1]
                if difference < params.DIFFERENCE_THRESH:
                    clump_counts += 1
                    if clump_counts == 1:
                        leading_counts += 1
                        if energies[i - 1] >= params.HIGH_ENERGY_LEAD_THRESH:  # Clump starts on the previous index
                            high_energy_lead += 1

                    # For the last count in the event
                    if i == self.stop - 1:
                        clump_counts += 1

                else:
                    # Adding to clumpiness when there's a clump of three or more
                    if clump_counts >= 3:
                        clumpiness += 1
                        clump_counts = 0

                    # Adding to the clumpiness when the gap between sufficient counts is greater than the threshold
                    if difference >= params.GAP_THRESH and clump_counts == 0:
                        clumpiness += 1

                    clump_counts = 0

        clumpiness /= self.length
        if leading_counts > 0:
            high_energy_lead /= leading_counts

        # Calculating the length subscore
        self.len_subscore = 1 / params.GOOD_LEN_THRESH * self.length

        # Calculating the clumpiness subscore
        self.clumpiness_subscore = 1 / (1 + np.e ** (40 * (clumpiness - params.CLUMPINESS_TOSSUP)))

        # Calculating high energy leading count subscore
        self.hel_subscore = np.e ** -(5.3 * high_energy_lead)

        # Getting weather subscore
        if detector.deployment['weather_station'] != '':
            local_date, local_time = tl.convert_to_local(detector.full_date_str, times[self.start])
            self.weather_subscore = tl.get_weather_conditions(detector, weather_cache, local_date, local_time)
        else:
            self.weather_subscore = -1

    def calculate_score(self, detector, weather_cache, times, energies):
        """Calculates the event's score (likelihood of being an interesting event) and updates its total_score
        attribute.

        Parameters
        ----------
        detector : tgfsearch.detectors.detector.Detector
            The Detector used to store all data and relevant information
        weather_cache : dict
            A cache containing weather info for the day being investigated. If no info has been added yet, this
            function will find and add it during the scoring process.
        times : numpy.ndarray
            An array containing times for each count.
        energies : np.array
            An array containing energies for each count.

        """

        assert (params.LEN_WEIGHT + params.CLUMPINESS_WEIGHT + params.HEL_WEIGHT + params.WEATHER_WEIGHT) == 1.0
        self._calculate_subscores(detector, weather_cache, times, energies)
        # If weather info couldn't be obtained, the weather subscore is removed and the remaining weights are adjusted
        # so that they stay proportional to one another
        if self.weather_subscore == -1:
            proportionality = 1 / (params.LEN_WEIGHT + params.CLUMPINESS_WEIGHT + params.HEL_WEIGHT)
            self.total_score = proportionality * (params.LEN_WEIGHT * self.len_subscore +
                                                  params.CLUMPINESS_WEIGHT * self.clumpiness_subscore +
                                                  params.HEL_WEIGHT * self.hel_subscore)
        else:
            self.total_score = (params.LEN_WEIGHT * self.len_subscore + 
                                params.CLUMPINESS_WEIGHT * self.clumpiness_subscore +
                                params.HEL_WEIGHT * self.hel_subscore +
                                params.WEATHER_WEIGHT * self.weather_subscore)
