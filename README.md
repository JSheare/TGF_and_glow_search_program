# TGF and Glow Search Program
## This program was written to search raw data from the instruments of the terrestrial gamma-ray flash group at UCSC.
## Currently supported detectors: THOR, GODOT, Santis Instrument (after adafruit update), Croatia instrument.

## **Requirements for use:**
- Python 3.10 (or later)
- All the following .py files:
  - search.py (core program)
  - search_module.py (program tools and functions)
  - search_classes.py (program objects)
  - searchGUI.py (program GUI)
  - DataReaderFinal.py (functions for reading and importing data files)
- All the following Python libraries:
  - scipy
  - numpy
  - pandas
  - matplotlib
  - psutil
  - Selenium
  - lxml

# **Instructions for Use:**
First, open a command line and navigate to the directory where all the .py files are located.

## **Using the GUI (Preferred):**
In the command line, type the following command:

    python3 searchGUI.py

Note that if you are running the program on Windows you may need to replace 'python3' with 'python'. This 
command will open the program's graphical user interface. 

To start a search, enter the first date in your search range in the entry box labeled 'Date one'. 

Next, enter the last date in your search range in the entry box labeled 'Date Two' (enter the same date as in the first 
box if you only want to search a single day). 

Afterward, specify the detector you want to search in the box labeled 'Detector' (a full list of supported detectors 
appears to the right of the box).

If you're running the program on Sol (UCSC TGF group computer) you're now ready to begin a search. Otherwise, 
you may need to specify an import directory for the data. This can be done by entering the data directory in the box
labeled 'Import Location' or by clicking on the 'Browse' button to the right of the box and navigating to the data
directory. 

Finally, to begin a search, simply click the 'Start' button. Output from the search will now appear onscreen, and
search results will be stored in a folder called 'Results' located inside the directory specified in the 
'Export Location' entry box (same as program location by default).

The search can be stopped at any time by clicking the 'Stop' button, and text can be cleared from onscreen by
clicking 'Clear Text'. You can also reset the entered dates/detector/modes all at once by clicking the 'Reset' button.

### **Program Modes:**
The program has several 'modes' that give it enhanced functionality. To run the program with any of the modes, simply
check the mode's corresponding tick box. Here's a list of each mode and their functions:
<br/>
#### **Normal Modes:**

- **'allscints' Mode** - This mode tells the program to run the short event search algorithm on all the scintillators 
(by default the program only searches the large plastic scintillator).


- **'template' Mode** - In order to accurately calibrate the large plastic scintillator, the program needs a template. New templates must be made
for each new detector location, and this mode allows you to make the template. To make templates, follow the instructions given by the program to adjust the vertical lines to match the locations of
the compton shoulders. A picture of what the template should look like is located in the program's GitHub repository:
https://github.com/JSheare/TGF_and_glow_search_program (look for 'lp calibration example.JPG').


- **'aircraft' Mode** - If the detector you are searching was deployed to an aircraft during the search period, this is 
the mode that you want to use. It adds a special check to the short event search algorithm that filters out false 
alarms caused by cosmic rays. It's recommended that the mode 'skcali' is used in conjunction with this.


- **'processed' Mode:** - In this mode (which is only available for GODOT) the program will search processed data instead 
of raw data. Note that this mode is largely deprecated and that processed data must exist for your specified search 
range for it to work properly.

#### **Developer Modes:**
The program also includes several bonus modes that are largely intended for developer use; 
they include:
- **'skcali' Mode** - Skips the scintillator calibration algorithm.


- **'skshort' Mode** - Skips the short event search algorithm.


- **'skglow' mode** - Skips the long event search algorithm.


- **'pickle' mode** - Serializes and saves the detector object for later use OR imports an 
  already-serialized detector object

## **Running the Program Through the Command Line:**

### **Single day:**
To search a single day only, enter a command of the following form:

    python3 search.py yymmdd yymmdd detector

where 'yy' corresponds to the year, 'mm' corresponds to the month, and 'dd' corresponds to the day.
<br/>
Note that if you are running the program on Windows you may need to replace 'python3' with 'python'.
<br/>
Here's an example for GODOT data on December 3rd, 2015:

    python3 search.py 151203 151203 GODOT

### **Date Range:**
To search a range of dates, follow the same instructions given for a single day, but replace the second 'yymmdd' 
with the last date in the desired range.
Here's an example for THOR5 data from July 1st, 2022 to August 31st, 2022:

    python3 search.py 220701 220831 THOR5

### **Program Modes:**
The program has several 'modes' that give it enhanced functionality. Here's a list of them all and how to use them:
<br/>
#### **'custom' Mode:**
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

    python3 search.py yymmdd yymmdd detector custom

#### **'allscints' Mode:**
This mode tells the program to run the short event search algorithm on all the scintillators 
(by default the program only searches the large plastic scintillator).
To run the program in this mode, run the program the same as above but add the word 'plastics' to the end:

    python3 search.py yymmdd yymmdd detector allscints

#### **'template' Mode:**
In order to accurately calibrate the large plastic scintillator, the program needs a template. New templates must be made
for each new detector location, and this mode allows you to make the template.
<br/>
To run the program in this mode, run the program the same as above but add the word 'template' to the end:

    python3 search.py yymmdd yymmdd detector template

To make templates, follow the instructions given by the program to adjust the vertical lines to match the locations of
the compton shoulders. A picture of what the template should look like is located in the program's GitHub repository:
https://github.com/JSheare/TGF_and_glow_search_program (look for 'lp calibration example.JPG').

#### **'aircraft' Mode:**
If the detector you are searching was deployed to an aircraft during the search period, this is the mode that you want 
to use. It adds a special check to the short event search algorithm that filters out false alarms caused by cosmic rays.
It's recommended that the mode 'skcali' is used in conjunction with this.
</br>
To run the program in this mode, run the program the same as above but add the word 'aircraft' to the end:

    python3 search.py yymmdd yymmdd detector aircraft

#### **'processed' Mode:**
In this mode (which is only available for GODOT) the program will search processed data instead of raw data. Note that
this mode is largely deprecated and that processed data must exist for your specified search range for it to work 
properly.
</br>
To run the program in this mode, run the program the same as above but add the word 'processed' to the end:

    python3 search.py yymmdd yymmdd detector processed

#### **Other Things to Know:**
It is possible to use as many of these modes in tandem as the user needs 
<br/>
(i.e. commands like this are possible):

    python3 search.py yymmdd yymmdd detector custom plastics template

### **Additional (developer) modes:**
The program also includes several bonus modes that are largely intended for developer use; 
they include:
- **'skcali' Mode** - Skips the scintillator calibration algorithm.


- **'skshort' Mode** - Skips the short event search algorithm.


- **'skglow' mode** - Skips the long event search algorithm.


- **'pickle' mode** - Serializes and saves the detector object for later use OR imports an 
  already-serialized detector object

### **Changing the export location of program results:**
In order to change the location where the program exports its results 
(i.e. short event scatter plots, long event histograms, scintillator spectra, logs)
follow these steps:
- Use a text editor or your favorite IDE to open search_module.py and look for the function 'results_loc'
- Once you've found the function, look for the variable called 'path' and enter your custom path inside the single 
  quotes to the right of the = sign
- The final result should look something like this (with example path /users/user/desktop):

    path = '/users/user/desktop'
- Save and close

The program should now save its results to the specified custom directory.
