# TGF and Glow Search Program
## This program was written to search raw data from the intstruments of the terrestrial gamma-ray flash group at UCSC.
## Currently supported detectors: THOR, GODOT, Santis Instrument (after adafruit update)
<br/>

##**Requirements for use:**
- python 3.6.6 (or later)
- All of the following .py files:
  - search.py (core program)
  - search_module.py (program tools and functions)
  - search_classes.py (program objects)
  - DataReaderFinal.py (functions for reading and importing data files)
- All of the following python libraries:
  - scipy
  - numpy
  - pandas
  - matplotlib
  - warnings
  - datetime
  - psutil

##**Instructions for Use:**
First, open a command line and navigate to the directory where search.py is located.

###**Single day:**
To search a single day only, enter a command of the following form:

    python search.py yymmdd yymmdd detector

where 'yy' corresponds to the year, 'mm' corresponds to the month, and 'dd' corresponds to the day.
<br/>
Note that if you have multiple python versions installed you may need to replace 'python' with 'python3'.
<br/>
Here's an example for GODOT data on December 3rd, 2015:

    python search.py 151203 151203 GODOT

###**Date Range:**
To search a range of dates, follow the same instructions given for a single day, but replace the second 'yymmdd' 
with the last date in the desired range.
Here's an example for THOR5 data from July 1st, 2022 to August 31st, 2022:

    python search.py 220701 220831 THOR5

###**Program Modes:**
The program has several 'modes' that give it enhanced functionality. Here's a list of them all and how to use them:
<br/>
####**'custom' Mode:**
This mode instructs the program to look for data in a custom, user-specified location (this feature is designed mainly 
for those using the program who are not on Sol, the primary host of the UCSC TGF groups' raw data).
<br/>
Before using this mode, the custom file path must first be specified. To do so, follow these steps:
- Use a text editor or your favorite IDE to open search_module.py and look for the function 'C_raw_data_loc'
- Once you've found the function, look for the variable called 'path' and enter your custom path inside the single 
  quotes to the right of the = sign
- The final result should look something like this (with example path /users/user/desktop):

    path = '/users/user/desktop'
- Remember to save before closing!

Once you've done this, run the program the same as above but add the word 'custom' to the end:

    python search.py yymmdd yymmdd detector custom

####**'plastics' Mode:**
This mode tells the program to run the short event search algorithm on all plastic scintillators 
(by default the program only searches the large plastic scintillator).
To run the program in this mode, run the program the same as above but add the word 'plastics' to the end:

    python search.py yymmdd yymmdd detector plastics

####**'template' Mode:**
In order to accurately calibrate the large plastic scintillator, the program needs a template. New templates must be made
for each new detector location, and this mode allows you to make the template.
<br/>
To run the program in this mode, run the program the same as above but add the word 'template' to the end:

    python search.py yymmdd yymmdd detector template

To make templates, follow the instructions given by the program to adjust the vertical lines to match the locations of
the compton shoulders. A picture of what the template should look like is located in the program's github repository:
https://github.com/JSheare/TGF_and_glow_search_program (look for 'lp calibration example.JPG').

####**Other Things to Know:**
It is possible to use as many of these modes in tandem as the user needs 
<br/>
(i.e. commands like this are possible):

    python search.py yymmdd yymmdd detector custom plastics template

###**Additional (developer) modes:**
The program also includes several bonus modes that can be used to skip various algorithms; 
they include:
- skcali: Skips the scintillator calibration algorithm.


- skshort: Skips the short event search algorithm.


- skglow: Skips the long event search algorithm.


- pickle: Pickles and saves the detector object for later use OR imports an 
  already-pickled detector object

###**Changing the export location of the program results:**
In order to change the location where the program exports its results 
(i.e. short event scatterplots, long event histograms, scintillator spectra, logs)
follow these steps:
- Use a text editor or your favorite IDE to open search_module.py and look for the function 'results_loc'
- Once you've found the function, look for the variable called 'path' and enter your custom path inside the single 
  quotes to the right of the = sign
- The final result should look something like this (with example path /users/user/desktop):

    path = '/users/user/desktop'
- Save and close

The program should now save its results to the custom directory specified.
