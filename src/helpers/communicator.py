"""A slightly modified version of the Event class from threading. Used by the GUI to communicate with the thread that
the search program is managed on"""
import threading as threading


class Communicator(threading.Event):
    def __init__(self):
        super().__init__()
        self.running = False
