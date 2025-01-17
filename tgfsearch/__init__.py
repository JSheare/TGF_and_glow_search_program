"""A package containing tools and scripts for searching raw data from the instruments of the TGF group at UCSC."""

# Everything here should be available at the top level of the package upon importing
from tgfsearch.tools import *
from tgfsearch.search import is_valid_detector, get_detector
from tgfsearch.DataReaderTimetrack2 import fileNameToData
from tgfsearch.detectors.croatia import Croatia
from tgfsearch.detectors.godot import Godot
from tgfsearch.detectors.santis import Santis
from tgfsearch.detectors.thor import Thor
from tgfsearch.detectors.adaptive_detector import AdaptiveDetector
from tgfsearch.events.shortevent import ShortEvent
from tgfsearch.events.longevent import LongEvent
from tgfsearch.helpers.reader import Reader
