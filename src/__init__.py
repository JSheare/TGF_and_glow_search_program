"""A package containing tools and scripts for searching raw data from the instruments of the TGF group at UCSC."""

# Everything here should be available at the top level of the package upon importing
from src.tools import *
from src.detectors.godot import Godot
from src.detectors.thor import Thor
from src.detectors.santis import Santis
from src.detectors.croatia import Croatia
from src.events.shortevent import ShortEvent
from src.events.longevent import LongEvent
from src.search import get_detector
from src.utilities.DataReaderTimetrack2 import fileNameToData
from src.utilities.DataReaderTimetrack2 import multifilesToData
