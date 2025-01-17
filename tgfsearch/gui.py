"""A graphical user interface for running the TGF search program."""
import multiprocessing as multiprocessing
import os as os
import platform as platform
import sys as sys
import threading as threading
import tkinter as tk
import traceback as traceback
from queue import Queue
from tkinter import filedialog
from tkinter import ttk

# Adds parent directory to sys.path. Necessary to make the imports below work when running this file as a script
parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import tgfsearch.tools as tl
from tgfsearch.search import is_valid_detector, mode_to_flag, program


# Redirects stdout and stderr from the search program
def search_program_wrapper(write, first_date, second_date, unit, mode_info):
    old_stdout_write = sys.stdout.write
    sys.stdout.write = write.send
    try:
        program(first_date, second_date, unit, mode_info)
    except Exception as ex:
        print('Search program terminated with the following error or warning:\n')
        # Removing the top layer of the traceback (which is just this function) and printing the remainder
        count = len(traceback.extract_tb(ex.__traceback__)) - 1
        print(traceback.format_exc(limit=-count))

    sys.stdout.write = old_stdout_write


# A helper class that keeps track of the required arguments for a single search
class SearchArgs:
    def __init__(self, first_date, second_date, detector, mode_info):
        self.first_date = first_date
        self.second_date = second_date
        self.detector = detector
        self.mode_info = mode_info

    def __str__(self):
        search_string = f'{self.first_date} {self.second_date} {self.detector}'
        for arg in self.mode_info:
            search_string += f' {arg}'

        return search_string

    def __hash__(self):
        return hash(tuple([self.first_date, self.second_date, self.detector] + self.mode_info))

    def __eq__(self, args2):
        return (self.first_date == args2.first_date and
                self.second_date == args2.second_date and
                self.detector == args2.detector and
                self.mode_info == args2.mode_info)


# A class for managing the search and keeping track of all the enqueued search information
class SearchManager:
    def __init__(self):
        self.mode_flags = {}
        self.search_queue = Queue()  # Queue that holds all the enqueued searches
        self.search_set = set()  # Set that keeps track of already-enqueued searches so that no duplicates are added
        self.lock = threading.Lock()  # Mutex to prevent race conditions with a search in progress
        self.event = threading.Event()  # Event for communicating with the search thread

    # Checks whether a search with the given dates and detector is a valid one
    @staticmethod
    def is_valid_search(first_date, second_date, detector, print_feedback=False):
        # Checks that both dates are digits in the proper format
        if not first_date.isdigit() or not second_date.isdigit() \
                or len(first_date) != 6 or len(second_date) != 6:
            if print_feedback:
                print('Error: not a valid date. BOTH dates must be in yymmdd format.')

            return False

        # Checks that both dates are sequential
        if int(first_date) > int(second_date):
            if print_feedback:
                print('Error: second date must be AFTER first date.')

            return False

        # Checks that a detector has been entered
        if detector == '':
            if print_feedback:
                print('Error: No detector specified.')

            return False
        elif not is_valid_detector(detector):
            if print_feedback:
                print('Error: Not a valid detector.')

            return False

        return True

    # Enables the given mode
    def add_mode(self, mode):
        if not self.lock.locked():
            self.mode_flags[mode] = True

    # Disables the given mode
    def remove_mode(self, mode):
        if not self.lock.locked() and mode in self.mode_flags:
            self.mode_flags[mode] = False

    # Enqueues a new search with the given parameters
    def enqueue(self, first_date, second_date, detector, import_loc, export_loc):
        if not self.lock.locked():
            if second_date == 'yymmdd' or second_date == '':
                second_date = first_date

            # If the search command is valid, sets up a SearchArgs object to store it
            if self.is_valid_search(first_date, second_date, detector, print_feedback=True):
                mode_info = []
                for mode in self.mode_flags:
                    if self.mode_flags[mode]:
                        mode_info.append(mode_to_flag(mode))

                mode_info.append('-c')
                if import_loc == '':
                    import_loc = 'none'

                if export_loc == '':
                    export_loc = 'none'

                mode_info.append(import_loc)
                mode_info.append(export_loc)
                search_args = SearchArgs(first_date, second_date, detector.upper(), mode_info)
                # Enqueues the search if it isn't a duplicate
                if search_args not in self.search_set:
                    # Reasoning behind 3: one for custom, the last two for custom import/export locations
                    modes_string = f' [{", ".join(mode_info[0:-3]).replace("-", "")}]' if len(mode_info) > 3 else ''
                    print(f'Enqueueing {tl.short_to_full_date(first_date)}'
                          f'{" - " + tl.short_to_full_date(second_date) if first_date != second_date else ""}'
                          f' on {detector.upper()}{modes_string}.')
                    self.search_queue.put(search_args)
                    self.search_set.add(search_args)

    # Runs all the enqueued searches
    def run(self):
        with self.lock:
            while not self.search_queue.empty():
                search_args = self.search_queue.get()
                self.search_set.remove(search_args)

                # Prints feedback about what date and modes were selected
                feedback_string = f'\nRunning search for {tl.short_to_full_date(search_args.first_date)}'
                if search_args.first_date != search_args.second_date:
                    feedback_string += f' - {tl.short_to_full_date(search_args.second_date)}'

                feedback_string += f' on {search_args.detector}.'
                print(feedback_string)
                # Reasoning behind 3: one for custom, the last two for custom import/export locations
                if len(search_args.mode_info) > 3:
                    print(f'This search will be run with the following modes: '
                          f'{", ".join(search_args.mode_info[0:-3]).replace("-", "")}.')

                # Runs the search program in a separate process and manages it
                read, write = multiprocessing.Pipe()
                process = multiprocessing.Process(target=search_program_wrapper,
                                                  args=(write, search_args.first_date, search_args.second_date,
                                                        search_args.detector, search_args.mode_info))
                process.start()
                while process.is_alive() and not self.event.is_set():
                    # Prints the processes' piped stdout
                    if read.poll():
                        print(read.recv(), end='')

                if process.is_alive():  # This will be executed if stop() is run
                    process.terminate()
                    break

                # In case there's still some strings left in the pipe
                while read.poll():
                    print(read.recv(), end='')

            self.event.clear()
            print('\nSearch Concluded.\n')

    # Stops the search if it's running
    def stop(self):
        if self.lock.locked():
            self.event.set()

    # Clears the search queue and resets the selected modes
    def reset(self):
        if not self.lock.locked():
            while not self.search_queue.empty():
                self.search_queue.get()

            self.search_set.clear()
            for mode in self.mode_flags:
                self.mode_flags[mode] = False


# A class implementing the search GUI window, all its widgets, and their associated functionality
class SearchWindow(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.master = master
        self.checkbox_variables = []

        # Making and placing the text box display
        self.text_box = tk.Text(self, height=30, width=100)
        self.text_box['state'] = tk.DISABLED
        self.text_box.grid(row=0, column=1, columnspan=3)

        self.text_box_label = tk.Label(self, text='Search Output')
        self.text_box_label.grid(row=1, column=2, pady=(5, 0))

        # Setting up the input frame
        self.input_frame = tk.Frame(self)
        self.input_frame.grid(row=2, column=1, rowspan=2)

        # Adding and placing the date/detector labels and entry boxes to the frame
        self.date_one_label = tk.Label(self.input_frame, text='Date One:')
        self.date_one_label.grid(row=0, column=0, pady=(5, 0))

        self.date_one_entry = tk.Entry(self.input_frame, width=15, borderwidth=5)
        self.date_one_entry.insert(0, 'yymmdd')
        self.date_one_entry.bind('<FocusIn>', lambda e: self._clear_ghost_text(self.date_one_entry, 'yymmdd'))
        self.date_one_entry.grid(row=1, column=0, pady=(5, 0))

        self.date_two_label = tk.Label(self.input_frame, text='Date Two:')
        self.date_two_label.grid(row=2, column=0, pady=(5, 0))

        self.date_two_entry = tk.Entry(self.input_frame, width=15, borderwidth=5)
        self.date_two_entry.insert(0, 'yymmdd')
        self.date_two_entry.bind('<FocusIn>', lambda e: self._clear_ghost_text(self.date_two_entry, 'yymmdd'))
        self.date_two_entry.grid(row=3, column=0, pady=(5, 0))

        self.detector_label = tk.Label(self.input_frame, text='Detector:')
        self.detector_label.grid(row=4, column=0, pady=(5, 0))

        self.detector_entry = tk.Entry(self.input_frame, width=15, borderwidth=5)
        self.detector_entry.grid(row=5, column=0, pady=(5, 0))

        # Setting up the search control frame
        self.search_frame = tk.Frame(self)
        self.search_frame.grid(row=2, column=2, rowspan=2)

        # Adding and placing the start, enqueue, and stop buttons, and the enqueue counter
        self.start_button = tk.Button(self.search_frame, height=3, width=20, text='Start', bg='white',
                                      command=self.start)
        self.start_button.grid(row=0, column=0, columnspan=2, pady=(5, 0))

        self.enqueue_button = tk.Button(self.search_frame, height=3, width=8, text='Enqueue', bg='white',
                                        command=self.enqueue)
        self.enqueue_button.grid(row=1, column=0, pady=(5, 0))

        self.stop_button = tk.Button(self.search_frame, height=3, width=8, text='Stop', bg='white',
                                     command=self.stop)
        self.stop_button.grid(row=1, column=1, pady=(5, 0))

        self.enqueue_label = tk.Label(self.search_frame, text='')
        self.enqueue_label.grid(row=2, column=0, columnspan=2, pady=(5, 0))

        # Setting up the display frame
        self.display_frame = tk.Frame(self)
        self.display_frame.grid(row=2, column=3)

        # Adding and placing the clear text and reset/ clear queue buttons
        self.clear_button = tk.Button(self.display_frame, height=3, width=8, text='Clear\nText', bg='white',
                                      command=self.clear)
        self.clear_button.grid(row=0, column=0, padx=(0, 4), pady=(5, 0))

        self.reset_button = tk.Button(self.display_frame, height=3, width=8, text='Reset/\nClear\nQueue', bg='white',
                                      command=self.reset)
        self.reset_button.grid(row=0, column=1, padx=(4, 0), pady=(5, 0))

        # Setting up the modes frame
        self.modes_frame = tk.Frame(self)
        self.modes_frame.grid(row=3, column=3, pady=(5, 0))

        # Adding and placing the regular mode label and checkboxes
        self.regular_cb_label = tk.Label(self.modes_frame, text='Modes:')
        self.regular_cb_label.grid(row=0, column=0, pady=(5, 0))

        ccb = tk.IntVar()
        self.combo_cb = tk.Checkbutton(self.modes_frame, text='combo', variable=ccb, onvalue=1, offvalue=0,
                                       command=lambda: self._check_uncheck(ccb, 'combo'))
        self.combo_cb.grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.checkbox_variables.append(ccb)

        ascb = tk.IntVar()
        self.allscints_cb = tk.Checkbutton(self.modes_frame, text='allscints', variable=ascb, onvalue=1, offvalue=0,
                                           command=lambda: self._check_uncheck(ascb, 'allscints'))
        self.allscints_cb.grid(row=2, column=0, sticky=tk.W, pady=(5, 0))
        self.checkbox_variables.append(ascb)

        acb = tk.IntVar()
        self.aircraft_cb = tk.Checkbutton(self.modes_frame, text='aircraft', variable=acb, onvalue=1, offvalue=0,
                                          command=lambda: self._check_uncheck(acb, 'aircraft'))
        self.aircraft_cb.grid(row=3, column=0, sticky=tk.W, pady=(5, 0))
        self.checkbox_variables.append(acb)

        # Adding and placing the developer mode label and checkboxes
        self.dev_cb_label = tk.Label(self.modes_frame, text='Dev Modes:')
        self.dev_cb_label.grid(row=0, column=1, pady=(5, 0))

        sscb = tk.IntVar()
        self.skshort_cb = tk.Checkbutton(self.modes_frame, text='skshort', variable=sscb, onvalue=1, offvalue=0,
                                         command=lambda: self._check_uncheck(sscb, 'skshort'))
        self.skshort_cb.grid(row=1, column=1, sticky=tk.W, pady=(5, 0))
        self.checkbox_variables.append(sscb)

        sgcb = tk.IntVar()
        self.skglow_cb = tk.Checkbutton(self.modes_frame, text='skglow', variable=sgcb, onvalue=1, offvalue=0,
                                        command=lambda: self._check_uncheck(sgcb, 'skglow'))
        self.skglow_cb.grid(row=2, column=1, sticky=tk.W, pady=(5, 0))
        self.checkbox_variables.append(sgcb)

        pcb = tk.IntVar()
        self.pickle_cb = tk.Checkbutton(self.modes_frame, text='pickle', variable=pcb, onvalue=1, offvalue=0,
                                        command=lambda: self._check_uncheck(pcb, 'pickle'))
        self.pickle_cb.grid(row=3, column=1, sticky=tk.W, pady=(5, 0))
        self.checkbox_variables.append(pcb)

        # Setting up the file import/export frame
        self.file_frame = tk.Frame(self)
        self.file_frame.grid(row=5, column=1, columnspan=3, pady=(10, 0))
        self.file_frame.columnconfigure(3, {'minsize': 30})

        # Adding and placing the custom import label, entry box, and file dialogue button
        self.import_label = tk.Label(self.file_frame, text='Import Location:')
        self.import_label.grid(row=0, column=0, columnspan=2, pady=(5, 0))

        self.import_entry = tk.Entry(self.file_frame, width=40, borderwidth=5)
        self.import_entry.grid(row=1, column=0, columnspan=2, pady=(5, 0))

        self.import_button = tk.Button(self.file_frame, width=6, height=2, text='Browse',
                                       command=lambda: self._select_dir(self.import_entry))
        self.import_button.grid(row=1, column=2, pady=(5, 0))

        # Adding and placing the custom export label, entry box, and file dialogue button
        self.export_label = tk.Label(self.file_frame, text='Export Location:')
        self.export_label.grid(row=0, column=4, columnspan=2, pady=(5, 0))

        self.export_entry = tk.Entry(self.file_frame, width=40, borderwidth=5)
        self.export_entry.grid(row=1, column=4, columnspan=2, pady=(5, 0))

        self.export_button = tk.Button(self.file_frame, width=6, height=2, text='Browse',
                                       command=lambda: self._select_dir(self.export_entry))
        self.export_button.grid(row=1, column=6, pady=(5, 0))

        if platform.system() == 'Windows':
            ttk.Separator(self, orient='horizontal').place(x=565, y=580, relwidth=0.25)  # Modes separator line
            ttk.Separator(self, orient='horizontal').place(x=0, y=703, relwidth=1.0)  # Import/export separator line
        else:
            ttk.Separator(self, orient='horizontal').place(x=593, y=625, relwidth=0.25)  # Modes separator line
            ttk.Separator(self, orient='horizontal').place(x=0, y=750, relwidth=1.0)  # Import/export separator line

        # Redirecting stdout to the GUI text box
        self.old_stdout_write = sys.stdout.write
        sys.stdout.write = self.write

        self.search_manager = SearchManager()
        self.search_thread = None

        # Starting the enqueue counter updater
        self.enqueued_counter_interval = 20  # interval at which search queue size is checked in milliseconds
        self._update_enqueued_counter()

    def __del__(self):
        sys.stdout.write = self.old_stdout_write  # Restoring stdout.write to what it was before

    # Creates a file dialogue and then puts the selected directory into the specified text entry box
    @staticmethod
    def _select_dir(entry_box):
        directory = filedialog.askdirectory(initialdir='/')
        entry_box.delete(0, 'end')
        entry_box.insert(0, directory)

    # Clears the sample text from the given text entry box
    @staticmethod
    def _clear_ghost_text(entry_box, ghost_text):
        current_text = entry_box.get()
        if current_text == ghost_text:
            entry_box.delete(0, 'end')

    # Updates the enqueued searches counter periodically
    def _update_enqueued_counter(self):
        self.enqueue_label['text'] = f'Searches\nEnqueued:\n{self.search_manager.search_queue.qsize()}'
        self.after(self.enqueued_counter_interval, self._update_enqueued_counter)

    # Adds/removes the given mode from the search arguments when the corresponding checkbox is checked/unchecked
    def _check_uncheck(self, var, mode):
        if var.get() == 1:
            self.search_manager.add_mode(mode)
        else:
            self.search_manager.remove_mode(mode)

    # Substitute function for sys.stdout.write that appends the given string to the GUI's big text box
    def write(self, input_str):
        self.text_box['state'] = tk.NORMAL
        self.text_box.insert('end', input_str, 'last_insert')
        self.text_box.yview(tk.END)
        self.text_box['state'] = tk.DISABLED

    # Enables/disables GUI widgets depending on the action parameter
    def _change_widgets(self, action):
        self.start_button['state'] = action
        self.enqueue_button['state'] = action
        self.date_one_entry['state'] = action
        self.date_two_entry['state'] = action
        self.detector_entry['state'] = action
        self.import_entry['state'] = action
        self.export_entry['state'] = action

        self.combo_cb['state'] = action
        self.allscints_cb['state'] = action
        self.aircraft_cb['state'] = action
        self.pickle_cb['state'] = action
        self.skshort_cb['state'] = action
        self.skglow_cb['state'] = action

    # Enables all checkboxes/buttons
    def enable_widgets(self):
        self._change_widgets(tk.NORMAL)

    # Disables all checkboxes/buttons
    def disable_widgets(self):
        self._change_widgets(tk.DISABLED)

    # Enqueues a new search based on the current contents of all the text entry boxes
    def enqueue(self):
        if self.search_thread is None:
            self.search_manager.enqueue(self.date_one_entry.get(), self.date_two_entry.get(), self.detector_entry.get(),
                                        self.import_entry.get(), self.export_entry.get())

    # Starts running the enqueued searches
    def start(self):
        self.enqueue()  # In case the current info hasn't been enqueued yet
        if self.search_thread is None and not self.search_manager.search_queue.empty():
            self.search_thread = threading.Thread(target=self._run, args=())
            self.search_thread.start()  # Running the search in another thread to prevent the GUI from locking up

    # Target function for search thread. Disables/Enables the GUI elements while the searches are running
    def _run(self):
        self.disable_widgets()
        self.search_manager.run()
        self.enable_widgets()
        self.search_thread = None

    # Stops the search if there is one
    def stop(self):
        if self.search_thread is not None:
            self.search_manager.stop()

    # Clears the big text box
    def clear(self):
        self.text_box['state'] = tk.NORMAL
        self.text_box.delete('1.0', 'end')
        self.text_box['state'] = tk.DISABLED

    # Stops the search and resets the GUI widgets back to their default states
    def reset(self):
        self.stop()
        self.search_manager.reset()  # Clearing the search queue and mode flags
        self.clear()

        self.date_one_entry.delete(0, 'end')
        self.date_one_entry.insert(0, 'yymmdd')
        self.date_two_entry.delete(0, 'end')
        self.date_two_entry.insert(0, 'yymmdd')
        self.detector_entry.delete(0, 'end')

        self.import_entry.delete(0, 'end')
        self.export_entry.delete(0, 'end')

        # Unchecking the checkboxes
        for variable in self.checkbox_variables:
            variable.set(0)


def main():
    # For running the program with pythonw (no terminal)
    if sys.stdout is None:
        sys.stdout = open(os.devnull, 'w')

    if sys.stderr is None:
        sys.stderr = open(os.devnull, 'w')

    root = tk.Tk()
    root.title('TGF Search')
    if platform.system() == 'Windows':
        root.geometry('1080x785')
    else:
        root.geometry('1080x845')

    gui = SearchWindow(root)
    gui.pack()
    root.mainloop()


if __name__ == '__main__':
    main()
