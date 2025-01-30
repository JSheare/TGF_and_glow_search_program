"""A helper class used to keep track of aligned traces"""


class TraceInfo:
    def __init__(self, trace_name, buff_no, times, energies):
        self.trace_name = trace_name
        self.buff_no = buff_no
        self.times = times
        self.energies = energies
