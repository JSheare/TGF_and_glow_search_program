### To-do list:
- [ ] Write duplicate file eliminator/sorter for general data
- [ ] Add full list of all locations ever for each detector
- [ ] Profile the program for bottlenecks
- [ ] Change detector object to accept yymmdd string instead of first second of day

### Completed
- [x] Fix channel-to-MeV converter function
- [x] Add timestamp to scatter plot title
- [x] Update scintillator calibration code for THOR
- [x] Find/fix beginning of day peak bug 
- [x] Replace all commands from calendar module with similar ones from datetime
- [x] Go through and do what all the improvement comments say
- [x] Review/change all for-loop variables 
- [x] Add docstrings to functions/classes
- [x] Figure out how to fix THOR5 MP eRC number change problem
- [x] Finish writing README.md
- [x] Add ability to request a range of dates to the image collector
- [x] Make event file generators more memory efficient
- [x] Add a dev mode that skips algorithm(s) of choice
- [x] Write an attribute updator method to detector class
- [x] Pickle detector object even when not in low memory mode?
- [x] Make a standard for adding more detectors
- [x] Make a child class of Detector called 'Chunk' for the low memory mode