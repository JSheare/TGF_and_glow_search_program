"""A package containing tools and scripts for searching raw data from the instruments of the TGF group at UCSC."""

# Everything here should be available at the top level of the package upon importing
from tgf_search.tools import *
from tgf_search.detectors.godot import Godot
from tgf_search.detectors.thor import Thor
from tgf_search.detectors.santis import Santis
from tgf_search.detectors.croatia import Croatia
from tgf_search.events.shortevent import ShortEvent
from tgf_search.events.longevent import LongEvent
from tgf_search.search import get_detector
