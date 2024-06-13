"""Parameters for use by the TGF search program and its modules."""

"""General Parameters"""
SEC_PER_DAY = 86400  # Number of seconds in a day
SEC_PER_HOUR = 3600  # Number of seconds in an hour
CENTURY = '20'  # The current century (numerically)
TWO_AM = 7200  # Number of seconds of the day corresponding to 2:00AM
OPERATING_MEMORY_ALLOWANCE = 100e6  # Maximum amount of memory that the program is allowed to use (in bytes)
# not including data
TOTAL_MEMORY_ALLOWANCE_FRAC = 0.25  # Fraction of available memory on the system that the program can use
MAX_CHUNKS = 16  # The maximum number of chunks allowed in low memory mode

K40_EDGE_ENERGY = 1.242  # The Compton edge energy for Potassium 40 (MeV)
T_EDGE_ENERGY = 2.381  # The Compton edge energy for Thorium (MeV)
K40_PHOTOPEAK_ENERGY = 1.46  # The Photo-peak energy for Potassium 40 (MeV)
T_PHOTOPEAK_ENERGY = 2.60  # The Photo-peak energy for Thorium (MeV)
T_K40_RATIO = 1.781  # The ratio of the decay energies of Thorium and Potassium 40. Used for NaI calibration.
NAI_CALI_RATIO_TOLERANCE = 0.001  # Used in NaI calibration. The max tolerance for the ratio between two peaks' energies
# to be considered "the same" as the desired ratio.


"""Trace-Related Parameters"""
LARGE_TRIGSPOT = 4092  # The usual trigspot for the large buffer
SMALL_TRIGSPOT = 1024  # The usual trigspot for the smaller buffers
TRIGGER_ABOVE_BASELINE = 2  # The number of energy channels above baseline needed to trigger a trace
ZERO_WEIGHT = 4  # The weight given pulse values of zero in the measurement of "bad" counts
GOOD_TRACE_THRESH = 1.5  # Traces with ratios of "good" and "bad" counts above this threshold are automatically passed
RISING_EDGE_MAX_SLOPE = 20  # The maximum average slope of a rising edge (which precedes saturation) in a trace
ROLLOVER_PERIOD = 2**36  # Maximum number of clock ticks before a rollover
T_STEP = 12.5e-9  # clock tick time resolution (80 MHz)
DT = 200e-9  # Time to extend sample on either side to let pulse shapes finish
TRACE_TRIGGER_THRESH = 5  # Threshold in mV for triggering a count from a trace
MV_PER_ADC = 0.2442002442002442  # Conversion for ADC scale (fixed)
DEADTIME_I = 24  # Holdoff (dead) time (units of samples)
INT_I = 24  # Integration time (units of samples)
PARTIAL_INT_I = 12  # Integration time (units of samples) for PSD partial integration
ENERGY_RESCALE = 1.0  # Possible energy rescaling for trace to counts
DEADTIME_EXTEND = 1  # Number of samples to extend deadtime by if not yet back to baseline after holdoff


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

# Scatter plot formatting parameters
SE_TIMESCALE_ONE = 5e-4  # Subplot 1 timescale, 500 Microseconds
SE_TIMESCALE_TWO = 0.005  # Subplot 2 timescale, 5 milliseconds
SE_TIMESCALE_THREE = 2  # Subplot 3 timescale, 2 seconds
NAI_COLOR = 'b'  # Color for NaI scintillator dots
SP_COLOR = 'm'  # Color for small plastic scintillator dots
MP_COLOR = 'g'  # Color for medium plastic scintillator dots
LP_COLOR = 'darkgoldenrod'  # Color for large plastic scintillator dots
DOT_ALPHA = 0.5  # Controls dot transparency


"""Long Event Search Parameters"""
# Search algorithm parameters
ENERGY_CUTOFF = 1.9  # MeV. All energies below this are cut out of the data during the long event search
BIN_SIZE = 4  # The size of each bin (in seconds)
FLAG_THRESH = 5  # The number of standard deviations above the mean at which a bin is flagged
LONG_EVENT_MIN_COUNTS = 1000  # Only in aircraft mode

# Rolling baseline parameters
WINDOW_SIZE = 20  # The number of bins in the window on each side
WINDOW_GAP = 5  # The number of bins between the center bin and the start of the window on each side

# Histogram subplot formatting parameters
LE_MAIN_BAR_COLOR = 'r'  # The color of the bars on the whole-day subplot
LE_MAIN_BAR_ALPHA = 0.5  # Controls the transparency of the bars in the whole-day subplot
LE_THRESH_LINE_COLOR = 'blue'  # Color of line representing triggering threshold
LE_SUBPLOT_PADDING = 20  # The max number of bins to pad a long event by on either side
LE_SUBPLOT_BAR_COLOR = 'c'  # Color of histogram bars in the subplots
LE_SUBPLOT_BAR_ALPHA = 0.5  # Controls the transparency of the bars in the subplots
