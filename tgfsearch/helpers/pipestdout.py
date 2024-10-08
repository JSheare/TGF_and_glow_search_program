"""A class for conveniently writing stdout to a multiprocessing connection object. Used by the GUI."""
import multiprocessing as multiprocessing
import sys as sys


class PipeStdout:
    def __init__(self, pipe):
        self.pipe = pipe
        self.stdout = sys.stdout
        sys.stdout = self

    def write(self, data):
        self.pipe.send(data)

    def flush(self):
        self.stdout.flush()

    def __del__(self):
        sys.stdout = self.stdout
