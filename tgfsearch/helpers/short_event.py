"""A helper class used to keep track of and process short events."""


class ShortEvent:
    def __init__(self, event_start, event_length, scintillator):
        self.start = int(event_start)
        self.length = int(event_length)
        self.stop = int(event_start + event_length)  # non-inclusive (in line with array slicing)
        self.scintillator = scintillator
        self.weather_conditions = 'no weather data'
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
        return f'ShortEvent({self.start}, {self.stop}, {self.scintillator})'

    # Debugging string dunder
    def __repr__(self):
        return self.__str__()
