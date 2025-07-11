"""A package containing tools and scripts for searching raw data from the instruments of the TGF group at UCSC."""

# Everything here should be available at the top level of the package upon importing
from tgfsearch.tools import *
from tgfsearch.search import is_valid_detector, get_detector
from tgfsearch.data_reader_annotated import fileNameToData
from tgfsearch.detectors.adaptive_detector import AdaptiveDetector
from tgfsearch.detectors.detector import Detector
from tgfsearch.helpers.reader import Reader
