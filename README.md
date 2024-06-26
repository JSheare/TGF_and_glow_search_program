# TGF and Glow Search Program
## A program that searches data from the instruments of the terrestrial gamma-ray flash group at UCSC.
## Currently supported detectors: Thor, Godot, Santis Instrument (after adafruit update), Croatia instrument.

## **Requirements for use:**
- Python 3.10 (or later)
- All the following Python libraries:
  - scipy
  - numpy
  - pandas
  - matplotlib
  - psutil
  - Selenium
  - lxml

# **Instructions for Use:**
To install the program, run the following pip command in the command line:

    pip install tgfsearch@git+https://github.com/JSheare/TGF_and_glow_search_program

## **Using the GUI (Preferred):**
In the command line, enter the following command:

    tgf-search

This will open the program's graphical user interface. 

To start a search, enter the first date in your search range in the entry box labeled 'Date One'. 

Next, enter the last date in your search range in the entry box labeled 'Date Two' (enter the same date as in the first 
box if you only want to search a single day). 

Afterward, specify the detector you want to search in the box labeled 'Detector' (a full list of supported detectors 
appears to the right of the box).

If you're running the program on Thrud (UCSC TGF group computer) you're now ready to begin a search. Otherwise, 
you may need to specify an import directory for the data. This can be done by entering the data directory in the box
labeled 'Import Location' or by clicking on the 'Browse' button to the right of the box and navigating to the data
directory. 

Finally, to begin a search, simply click the 'Start' button. Output from the search will now appear onscreen, and
search results will be stored in a folder called 'Results' located inside the directory specified in the 
'Export Location' entry box (present working directory by default).

The search can be stopped at any time by clicking the 'Stop' button, and text can be cleared from onscreen by
clicking 'Clear Text'. You can also reset the entered dates/detector/modes all at once by clicking the 'Reset' button.

### **Program Modes:**
The program has several 'modes' that give it enhanced functionality. To run the program with any of the modes, simply
check the mode's corresponding tick box. Here's a list of each mode and its function:
<br/>
#### **Normal Modes:**

- **'allscints' Mode** - This mode tells the program to run the short event search algorithm on all the scintillators 
individually (by default the program only searches the large plastic scintillator).


- **'template' Mode** - In order to accurately calibrate the large plastic scintillator, the program needs a template. 
New templates must be made for each new detector location, and this mode allows you to make the template. To make 
templates, follow the instructions given by the program and adjust the vertical lines to match the locations of
the compton shoulders. A picture of what the template should look like is located in the program's GitHub repository:
https://github.com/JSheare/TGF_and_glow_search_program (look for 'lp calibration example.JPG').


- **'aircraft' Mode** - If the detector you are searching was deployed to an aircraft during the search period, this is 
the mode that you want to use. It adds a special check to the short event search algorithm that filters out false 
alarms caused by cosmic rays. It's recommended that the mode 'skcali' is used in conjunction with this.


- **'combo' Mode:** - This mode instructs the program to combine data from all available scintillators 
when looking for TGFs.

#### **Developer Modes:**
The program also includes several bonus modes that are largely intended for developer use; 
they include:
- **'skcali' Mode** - Skips the scintillator calibration algorithm.


- **'skshort' Mode** - Skips the short event search algorithm.


- **'skglow' mode** - Skips the long event search algorithm.


- **'pickle' mode** - Serializes and saves imported data for later use OR imports 
  previously-serialized data.

### **Enqueuing Multiple Searches:**
If desired, several searches can be enqueued and then executed sequentially.
To do this, simply fill out the search fields as described above and then press the 'Enqueue' button
to the left of the display screen. Once you've added all of your desired searches to the queue,
click 'Start' and the program will execute them all in the order in which they were entered.

To clear the queue without running any searches, click the 'Reset' button.

## **Running the Program Through the Command Line:**

### **Single day:**
To search a single day only, enter a command of the following form:

    tgf-search-cl yymmdd yymmdd detector

where 'yy' corresponds to the year, 'mm' corresponds to the month, and 'dd' corresponds to the day.
<br/>
Here's an example for GODOT data on December 3rd, 2015:

    tgf-search-cl 151203 151203 GODOT

### **Date Range:**
To search a range of dates, follow the same instructions given for a single day, but replace the second 'yymmdd' 
with the last date in the desired range.
Here's an example for THOR5 data from July 1st, 2022 to August 31st, 2022:

    tgf-search-cl 220701 220831 THOR5

### **Program Modes:**
The program has several 'modes' that give it enhanced functionality. Here's a list of them all and how to use them:
<br/>
#### **'combo' Mode:**
This mode instructs the program to combine data from all available scintillators when looking for TGFs.
To use the program in this mode, enter a command the same way as above but add the flag '--combo' to the end:

    tgf-search-cl yymmdd yymmdd detector --combo

#### **'allscints' Mode:**
This mode tells the program to run the short event search algorithm on all the scintillators 
(by default the program only searches the large plastic scintillator).
To run the program in this mode, run the program the same as above but add the flag '--plastics' to the end:

    tgf-search-cl yymmdd yymmdd detector --allscints

#### **'template' Mode:**
In order to accurately calibrate the large plastic scintillator, the program needs a template. 
New templates must be made for each new detector location, and this mode allows you to make the template.
<br/>
To run the program in this mode, run the program the same as above but add the flag '--template' to the end:

    tgf-search-cl yymmdd yymmdd detector --template

To make templates, follow the instructions given by the program to adjust the vertical lines to match the locations of
the compton shoulders. A picture of what the template should look like is located in the program's GitHub repository:
https://github.com/JSheare/TGF_and_glow_search_program (look for 'lp calibration example.JPG').

#### **'aircraft' Mode:**
If the detector you are searching was deployed to an aircraft during the search period, this is the mode that you want 
to use. It adds a special check to the short event search algorithm that filters out false alarms caused by cosmic rays.
It's recommended that the mode 'skcali' is used in conjunction with this.
</br>
To run the program in this mode, run the program the same as above but add the flag '--aircraft' to the end:

    tgf-search-cl yymmdd yymmdd detector --aircraft

#### **'processed' Mode:**
In this mode (which is only available for GODOT) the program will search processed data instead of raw data. Note that
this mode is largely deprecated and that processed data must exist for your specified search range for it to work 
properly.
</br>
To run the program in this mode, run the program the same as above but add the flag '-p' to the end:

    tgf-search-cl yymmdd yymmdd detector -p

#### **'custom' Mode:**
This mode can be used to specify custom data import and result export locations. To use it, enter a command
of the following form:

    tgf-search-cl yymmdd yymmdd detector -c import_directory export_directory

Note: if you don't wish to specify one of the directories, simply use the word 'none' instead. Here's an example
where we omit a custom import directory:

    tgf-search-cl yymmdd yymmdd detector -c none export_directory

#### **Other Things to Know:**
It is possible to use as many of these modes in tandem as the user needs 
<br/>
(i.e. commands like this are possible):

    tgf-search-cl yymmdd yymmdd detector -c --allscints --template

### **Additional (developer) modes:**
The program also includes several bonus modes that are largely intended for developer use; 
they include:
- **'skcali' Mode** - Skips the scintillator calibration algorithm. Flag: --skcali


- **'skshort' Mode** - Skips the short event search algorithm. Flag: --skshort


- **'skglow' mode** - Skips the long event search algorithm. Flag: --skglow


- **'pickle' mode** - Serializes and saves imported data for later use OR imports 
  previously-serialized data. Flag: --pickle

# **Using Data Tools:**
In addition to the search program, the package also includes data handling tools that you can use to build your own 
programs and scripts. This section will serve as an example for how to use many of these tools, with special emphasis 
placed on Detector objects, as they are the main focus of the package. 

To get started, import the package:
```python3
import numpy as np  # for parts of the example
import tgfsearch as tgf
```
Just to make sure we've covered all our bases, for those who are already familiar with the single-file data reader, it's
included in the package and can be used like so:
```python3
passtime = {'lastsod': -1.0, 'ppssod': -1.0, 'lastunix': -1.0, 'ppsunix': -1.0, 'lastwc': 0,
            'ppswc': 0, 'hz': 8e7, 'started': 0}
data = tgf.fileNameToData('/media/tgfdata/Detectors/THOR/THOR1/220831/eRC4195lpl_lm_220831_235435.txt.gz', passtime)
# where data is just a dataframe like you're probably used to.
```

### **Basic Detector Class Usage:**
Detector is a class that's designed to quickly and efficiently import and provide access to all data (list mode and 
trace) for a single day. Throughout these examples, we'll be using August 31st, 2022 on Thor 1 as a test day.

The best way to get a detector object of the correct type is to use the provided function get_detector():
```python3
detector = tgf.get_detector('THOR1', '220831')
```
If you're using the package on Thrud, you're already good to go. But in case you aren't, or you're trying to import data
from somewhere other than the main data drive, you can also set a custom import location:
```python3
detector.set_import_loc('/home/user/THOR1/220831')
```
Additionally, some of Detector's member functions are capable of exporting files for various purposes. By default, these
files are exported to your present working directory, but you can also specify a specific directory like so:
```python3
detector.set_results_loc('/home/user')
```
Now we're ready to import data, and we do so like this:
```python3
detector.import_data()
```
Note that this process can take several minutes depending on the size of the data set. By default, both list mode data 
and trace data are imported, but this can be changed using the two optional boolean parameters import_lm and 
import_traces.

The Detector will attempt to import any data it finds in its specified import directory. But, in case the data is 
missing or couldn't be imported properly, you can also check to see if it's present:
```python3
for scintillator in detector:
    if not detector.data_present_in(scintillator):  # Checking that list mode data is present (default)
        print(f'List mode data missing for {scintillator}.')

    if not detector.data_present_int(scintillator, data_type='trace'):  # Checking that trace data is present
        print(f'Trace data missing for {scintillator}.')
```

### **Getting and Setting Basic Information for Detector:**
The Detector class has several useful class attributes containing basic information that are worth knowing about:
```python3
unit = detector.unit  # The name of the detector (string).
date_str = detector.date_str  # The date in yymmdd format (string).
full_date_Str = detector.full_date_str  # The date in yyyy-mm-dd format (string).
first_sec = detector.first_sec  # The first second of the day in epoch time.
scint_list = detector.scint_list  # A list of abbreviations corresponding to the instrument's scintillators.
```
These are the most likely to be useful, but there are several others too. See the Detector documentation (linked at the 
bottom) for a full list.

If the information you're looking for is particular to a single scintillator, you can use the get_attribute() method:
```python3
eRC = detector.get_attribute('LP', 'eRC')  # Getting the serial number of the large plastic scintillator.
lm_filelist = detector.get_attribute('LP', 'lm_filelist')  # Getting a list of all large plastic list mode files.
trace_filelist = detector.get_attribute('LP', 'trace_filelist')  # Getting a list of all large plastic trace files.
lm_data = detector.get_attribute('LP', 'lm_frame')  # Getting all the large plastic list mode data for the whole day.
```
As above, these are the most common attributes, but you can request any on the following list:
- name - The scintillator's name (abbreviated).
- eRC - The scintillator's serial number.
- lm_frame - A pandas dataframe containing all the scintillator's list mode data.
- calibration_bins - A list containing the energy bins corresponding to Compton edges/photo peaks used for calibration.
- calibration_energies - A list containing the energies of Compton edges/photo peaks used for calibration.
- lm_files - A list of list mode files for the day.
- lm_file_ranges - A list of lists. Each sublist contains a pair of numbers corresponding to the first and last second
        in each list mode file.
- lm_file_indices - A dictionary of lists. Each list contains the indices needed to slice data for a particular file
        out of lm_frame.
- trace_filelist - A list of trace files for the day.
- traces - A dictionary containing trace data for each of the day's trace files.
- passtime - A dictionary containing information needed to import the subsequent list mode files properly
        (if applicable).

Additionally, you can also assign new information to these attributes using the Detector's set_attribute() method:
```python3
# Setting the filelist for the large plastic scintillator to a single file.
detector.set_attribute('LP', 'lm_filelist',
                       ['/media/tgfdata/Detectors/THOR/THOR1/220831/eRC4195lpl_lm_220831_235435.txt.gz'])
```
Note that the new information must be of the same type as the old information or an error will be raised.

### **Getting and Setting List Mode Data Columns:**
Detector has a built-in method for getting individual columns of list mode data as numpy arrays.
```python3
time = detector.get_lm_data('LP', 'SecondsOfDay')  # Getting the large plastic SecondsOfDay column as a numpy array.
energy = detector.get_lm_data('LP', 'energy')  # Getting the large plastic energy column as a numpy array.
wc = detector.get_lm_data('LP', 'wc')  # Getting the large plastic wc (wallclock) column as a numpy array.
```
You can also request a column from only a single data file like so:
```python3
file_time = detector.get_lm_data('LP', 'SecondsOfDay', file_name='/media/tgfdata/Detectors/THOR/THOR1/220831/'
                                                                 'eRC4195lpl_lm_220831_235435.txt.gz')
```
If you're interested in the energy column, you can also optionally convert the energies to units of MeV. Note that the 
Detector must first be calibrated for this to work, though.
```python3
# Additional comment: only available for LP and NaI at the moment.

detector.calibrate()  # calibrating the Detector
mev_energy = detector.get_lm_data('LP', 'energy', to_mev=True)
```
Finally, you can also set the data for a particular column using the set_lm_data() method:
```python3
detector.set_lm_data('LP', 'energy', np.array([1, 2, 3, 4, 5, 6, 7]))  # Example array
```
Note that the new data must be of the same length as the old data or an error will be raised.

### **Finding Matching List Mode Files and Trace Files:**
Detector contains two main member functions for matching count times with their corresponding list mode files or trace 
files. Say, for example, we have a large plastic count at 73000 seconds of day that we want to match:
```python3
# Matching the count to the name (as a string) of the list mode file it came from.
matching_lm_file = detector.find_lm_file('LP', 73000)

# Matching the count to a list of trace file names (strings) that it *could* have come from.
matching_traces = detector.find_matching_traces('LP', 73000)
```

### **Getting Specific Files' Data:**
Detector has two member functions for getting all data from specific files: one is for list mode files, the other is for
trace files. All you need to provide is the name of the scintillator the file is from and the name of the file itself:
```python3
lm_file_data = detector.get_lm_file('LP', '/media/tgfdata/Detectors/THOR/THOR1/220831/'
                                          'eRC4195lpl_lm_220831_235435.txt.gz')

trace_file_data = detector.get_trace('LP', '/media/tgfdata/Detectors/THOR/THOR1/220831/'
                                           'eRC4195lpl_xtr_220831_191800.xtr.gz')
# These functions return dataframes.
```

### **Pickling and Unpickling Detectors:**
A Detector object can be pickled (saved to a file, data and all) and unpickled using the following functions:
```python3
# Pickling our detector object in a file named "my_detector.pickle" to a results folder in the directory "/home/user"
# Note that path is an optional variable, and if left unspecified the detector will be pickled to its results directory.
pickle_path = tgf.pickle_detector(detector, 'my_detector', path='/home/user')

detector = tgf.unpickle(pickle_path)
```
These functions can be useful if you're using the same data set over and over again and want to avoid the lengthy
import procedure.

### **Splicing Two Detectors Together:**
Data from two different Detectors can be combined using the splice() method, yielding a new Detector with
a combination of both datasets:
```python3
# Making another Detector for the purpose of the example
detector2 = tgf.get_detector('THOR1', '220830')
detector2.set_import_loc('/home/user/THOR1/220830')
detector2.import_data()

# Calling the splice() method returns a new Detector
detector3 = detector1.splice(detector2)

# If you want, you can also use the normal addition syntax to do a splice
# detector3 = detector + detector2

# or += to overwrite an existing detector with the new, spliced one
# detector += detector2
```
The new Detector (detector3) will contain the seamlessly-combined data from both detector1 and detector2.
There are a few things to note about the new Detector:
- It will use the *earlier* of the two Detectors' dates as its own. For example, because detector has the date
      220831, and detector2 has the date 220830, detector3 will use the earlier 220830.
- Time data in the new Detector will be updated to reflect the difference in date between the two provided
      Detectors. For example, detector and detector2 are separated by a day, so all the data from the later
      Detector (detector in this case) will be corrected one day forward.
- Trying to import data with the new Detector won't work and will instead raise an error.

If you ever need to get a list of dates stored in a Detector, you can do so with Detector's dates_stored attribute:
```python3
dates_stored = detector3.dates_stored  # In this example, the list would be ['220830', '220831']
```

### **Trace Filtering and Alignment:**
The package contains several functions for filtering and aligning traces that you may find useful.

```python3
# For single traces. is_good will be True if the trace is good.
is_good = tgf.is_good_trace(trace_file_data)

# For a group of traces stored in a Detector. The result will be a list containing the names of traces that are
# likely to be interesting for the requested scintillator (the large plastic in this example).
good_trace_names = tgf.filter_traces(detector, 'LP')

if is_good:
    # Aligning the trace with its corresponding list mode data. This will return two numpy arrays: one with the
    # correctly-aligned times and another with the magnitude-corrected trace energies.
    trace_times, trace_energies = tgf.align_trace(trace_file_data, lm_file_data)
```

### **Clearing Data from Detector:**
Detector also has a method for clearing all currently stored data:
```python3
detector.clear()
```
Additionally, you can also opt to clear all stored data BUT leave file lists alone:
```python3
detector.clear(clear_filelists=False)
```

### **Detector Methods and Deep Copying:**
To make sure that data can't be modified by accident, several of Detector's getter and setter functions return/make deep
copies of the resource being requested/assigned. This behavior can be overridden using the optional boolean parameter 
'deepcopy' for the following member functions:
- get_attribute()
- set_attribute()
- get_lm_file()
- get_trace()

For example, the following code will return the actual list mode file dataframe stored in Detector rather than a copy:
```python3
actual_lm_data = detector.get_lm_file('LP', '/media/tgfdata/Detectors/THOR/THOR1/220831/'
                                            'eRC4195lpl_lm_220831_235435.txt.gz', deepcopy=False)
```
This can be performant in some cases, but be careful not to modify your data by accident!

## **List of Common Functions and Methods:**

### **Detector:**
- Detector.import_data() <br/> 
  - Imports data from data files into arrays and then updates them into the detector's scintillator objects.
- Detector.get_import_loc() <br/>
  - Returns the directory where data will be imported from.
- Detector.set_import_loc()
  - Sets the directory where data will be imported from.
- Detector.get_results_loc()
  - Returns the directory where all results will be stored.
- Detector.set_results_loc()
  - Sets the directory where all results will be stored.
- Detector.is_named()
  - Returns True if the Detector has the same name as the passed string.
- Detector.data_present_in()
  - Returns True if data is present for the specified scintillator and False otherwise.
- Detector.get_attribute()
  - Returns the requested attribute for a particular scintillator.
- Detector.set_attribute()
  - Updates the requested attribute for a particular scintillator.
- Detector.get_lm_data()
  - Returns a single column of list mode data as a numpy array.
- Detector.set_lm_data()
  - Sets a single column of list mode data to the new data specified.
- Detector.find_lm_file()
  - Returns the name of the list mode file that the given count occurred in.
- Detector.get_lm_file()
  - Returns the list mode data for the specified list mode file.
- Detector.get_trace()
  - Returns the trace data for the given scintillator and trace name.
- Detector.get_trace_names()
  - Returns a list of names of the traces that are currently being stored.
- Detector.find_matching_traces()
  - Finds the traces that could be a match for the given count (if they exist).
- Detector.get_fileset_size()
  - Returns the total size (in bytes) of all the currently listed files for the day.
- Detector.clear()
  - Clears all data currently stored in the Detector.
- Detector.calibrate()
  - Makes energy spectra histograms and calibrates the large plastic and sodium iodide scintillators.

Full Detector documentation for these methods (and more) can be found [here](https://html-preview.github.io/?url=https://github.com/JSheare/TGF_and_glow_search_program/blob/master/docs/detector.html).

### **Tools:**
- tgf.print_logger()
  - Prints the specified string to both stdout and the specified file.
- tgf.days_per_month()
  - Returns the number of days in the requested month based on the year.
- tgf.roll_date_forward()
  - Returns the calendar date after the one given as an argument.
- tgf.roll_date_backward()
  - Returns the calendar date before the one given as an argument.
- tgf.make_date_list()
  - Makes a list of dates from first_date to second_date (inclusive).
- tgf.full_date_to_short()
  - Converts a date string of the form yyyy-mm-dd to the form yymmdd.
- tgf.short_to_full_date()
  - Converts a date string of the form yymmdd to the form yyyy-mm-dd.
- tgf.get_first_sec()
  - Converts the given date string (in yymmdd format) to its first second in EPOCH time.
- tgf.pickle_detector()
  - Pickles Detectors.
- tgf.unpickle_detector()
  - Unpickles Detectors.
- tgf.filter_files()
  - Returns an ordered list of files with duplicate/invalid files filtered out.
- tgf.separate_files()
  - Returns a pair of ordered lists: one with list mode files, the other with trace files.
- tgf.scrape_weather()
  - Scrapes weather from weather underground and returns the results as a pandas data frame.
- tgf.combine_data()
  - Combines data from all scintillators into one set of arrays.
- tgf.separate_data()
  - Separates combined data from multiple scintillators into separate data for each scintillator.
- tgf.is_good_trace()
  - Returns True if the given trace is likely to be interesting and False otherwise.
- tgf.filter_traces()
  - Returns a list of traces that are likely to be interesting for the given scintillator.
- tgf.align_trace()
  - Aligns the given trace with the given list mode data.

Full tools documentation for these functions (and others) can be found [here](https://html-preview.github.io/?url=https://github.com/JSheare/TGF_and_glow_search_program/blob/master/docs/tools.html).

## **Other Documentation:**
- [Scintillator](https://html-preview.github.io/?url=https://github.com/JSheare/TGF_and_glow_search_program/blob/master/docs/scintillator.html)
- [Short Event Helper](https://html-preview.github.io/?url=https://github.com/JSheare/TGF_and_glow_search_program/blob/master/docs/shortevent.html)
- [Long Event Helper](https://html-preview.github.io/?url=https://github.com/JSheare/TGF_and_glow_search_program/blob/master/docs/longevent.html)
