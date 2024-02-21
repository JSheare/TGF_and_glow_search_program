"""Parameters for use by the TGF search program and it's modules."""

"""General Parameters"""
SEC_PER_DAY = 86400  # Number of seconds in a day
SEC_PER_HOUR = 3600  # Number of seconds in an hour
CENTURY = '20'  # The current century (numerically)
TWO_AM = 7200  # Number of seconds of the day corresponding to 2:00AM
OPERATING_MEMORY_ALLOWANCE = 100e6  # Maximum amount of memory that the program is allowed to use (in bytes)
# not including data
TOTAL_MEMORY_ALLOWANCE_FRAC = 0.25  # Fraction of available memory on the system that the program can use
MAX_CHUNKS = 16  # The maximum number of chunks allowed in low memory mode


"""Short Event Search Parameters"""
# Search algorithm parameters
NORMAL_ROLLGAP = 4  # The normal short event search rollgap
COMBO_ROLLGAP = 12  # the rollgap used in combo mode
AIRCRAFT_ROLLGAP = 18  # the rollgap used in aircraft mode
SHORT_EVENT_TIME_SPACING = 1e-3  # 1 millisecond
SHORT_EVENT_MIN_COUNTS = 10  # The minimum number of counts that a short event needs to be
MAX_PLOTS_PER_SCINT = 1000  # The maximum number of scatter plots/event files that can be made per scintillator

# Filter/ranking parameters
GOOD_LEN_THRESH = 30  # The number of counts at which a short event becomes *definitely* interesting

#   Noise filter parameters
CHANNEL_RANGE_WIDTH = 300  # The width of channels examined by the low/high energy ratio filter
CHANNEL_SEPARATION = 100  # The number of channels separating the low/high channel ranges
LOW_CHANNEL_START = 200  # The starting channel of the low energy channel range
CHANNEL_RATIO = 0.5  # The low/high energy channel range ratio required by the low/high energy ratio filter
MIN_NOISE_COUNTS = 3  # The minimum number of non-noise counts required for an event
NOISE_CUTOFF_ENERGY = 300  # The threshold for what's considered a noise/non-noise count

#   Successive CRS filter/clumpiness parameters
DIFFERENCE_THRESH = 2e-6  # The maximum time separation between two counts in a single clump
GAP_THRESH = 10e-6  # The minimum time separation between two counts of two different clumps
CLUMPINESS_THRESH = 0.27  # The clumpiness above which an event is probably just a successive CRS
CLUMPINESS_TOSSUP = 0.2  # The clumpiness at which an event could either be successive CRS or a real event

#   High energy lead parameters
HIGH_ENERGY_LEAD_THRESH = 15000  # The cutoff for what is considered a high-energy count

#   Ranking weight parameters (should add up to 1)
LEN_WEIGHT = 0.3
CLUMPINESS_WEIGHT = 0.2
HEL_WEIGHT = 0.2
WEATHER_WEIGHT = 0.3


"""Long Event Search Parameters"""
# Search algorithm parameters
ENERGY_CUTOFF = 1.9  # MeV. All energies below this are cut out of the data during the long event search
BIN_SIZE = 4  # The size of each bin (in seconds)
FLAG_THRESH = 5  # The number of standard deviations above the mean at which a bin is flagged
LONG_EVENT_MIN_COUNTS = 1000  # Only in aircraft mode

# Rolling baseline parameters
WINDOW_SIZE = 20  # The number of bins in the window on each side
WINDOW_GAP = 5  # The number of bins between the center bin and the start of the window on each side
