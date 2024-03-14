"""A graphical user interface for running the TGF search program."""
import tkinter as tk
import sys as sys
import os as os
import threading as threading
import multiprocessing as multiprocessing
import traceback
from queue import Queue
from tkinter import filedialog

# Adds parent directory to sys.path. Necessary to make the imports below work when running this file as a script
parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import tgfsearch.tools as tl
from tgfsearch.search import program


# A slightly modified version of the Event class from threading.
# Used to communicate with the thread that the search program is managed on
class Communicator(threading.Event):
    def __init__(self):
        super().__init__()
        self.running = False


# A class for conveniently writing stdout to a multiprocessing connection object
class PipeStdout:
    def __init__(self, pipe):
        self.pipe = pipe
        self.stdout = sys.stdout
        sys.stdout = self

    def write(self, data):
        self.pipe.send(data)

    def flush(self):
        self.stdout.flush()

    def __del__(self):
        sys.stdout = self.stdout


# Creates a directory selection dialogue box and then puts the selected directory in the specified text entry box
def select_dir(entry_box):
    directory = filedialog.askdirectory(initialdir='/')
    entry_box.delete(0, 'end')
    entry_box.insert(0, directory)


# Clears the sample text from the date entry boxes when they are clicked
def clear_ghost_text(entrybox, ghost_text):
    current_text = entrybox.get()
    if current_text == ghost_text:
        entrybox.delete(0, 'end')


# Appends/removes mode from mode list when the corresponding tick box is ticked/unticked
def tick_untick(var, modes, mode):
    if var.get() == 1:
        modes.append(mode)
    else:
        modes.remove(mode)


# Clears any text from the big text box
def clear_box(gui):
    text_box = gui.nametowidget('text_box')
    text_box['state'] = tk.NORMAL
    text_box.delete('1.0', 'end')
    text_box['state'] = tk.DISABLED


# Enables/disables gui elements depending on the action parameter
def change_elements(gui, action):
    gui.nametowidget('start_button')['state'] = action
    gui.nametowidget('enqueue_button')['state'] = action
    gui.nametowidget('date_one')['state'] = action
    gui.nametowidget('date_two')['state'] = action
    gui.nametowidget('detector_entrybox')['state'] = action
    gui.nametowidget('custom_entrybox')['state'] = action
    gui.nametowidget('results_entrybox')['state'] = action

    gui.nametowidget('allscints')['state'] = action
    gui.nametowidget('template')['state'] = action
    gui.nametowidget('aircraft')['state'] = action
    gui.nametowidget('combo')['state'] = action
    gui.nametowidget('skcali')['state'] = action
    gui.nametowidget('skshort')['state'] = action
    gui.nametowidget('skglow')['state'] = action
    gui.nametowidget('pickle')['state'] = action


# Disables all checkboxes/buttons
def disable_elements(gui):
    change_elements(gui, tk.DISABLED)


# Enables all checkboxes/buttons
def enable_elements(gui):
    change_elements(gui, tk.NORMAL)


# Checks whether a search command is valid
def is_valid_search(first_date, second_date, detector):
    # Checks that both dates are digits in the proper format
    if not first_date.isdigit() or not second_date.isdigit() \
            or len(first_date) != 6 or len(second_date) != 6:
        print('Error: not a valid date. Please enter BOTH dates in yymmdd format.')
        return False

    # Checks that both dates are sequential
    if int(first_date) > int(second_date):
        print('Error: second date must be AFTER first date.')
        return False

    # Checks that a detector has been entered
    if detector == '':
        print('Error: Please enter a detector.')
        return False

    return True


# Checks whether a queue contains the specified item
def queue_contains(item, queue):
    queue_len = queue.qsize()
    in_queue = False
    for i in range(queue_len):
        temp_item = queue.get_nowait()
        if not in_queue:
            if temp_item == item:
                in_queue = True

        queue.put_nowait(temp_item)

    return in_queue


def enqueue(gui, search_queue, program_modes):
    first_date = gui.nametowidget('date_one').get()
    second_date = gui.nametowidget('date_two').get()
    detector = gui.nametowidget('detector_entrybox').get()
    # If the search command is valid, constructs the command and adds it to the queue
    if is_valid_search(first_date, second_date, detector):
        command = [first_date, second_date, detector.upper(), 'gui']
        for mode in program_modes:
            command.append(mode)

        custom_results_dir = gui.nametowidget('results_entrybox').get() if (
                gui.nametowidget('results_entrybox').get() != '') else 'none'
        custom_import_dir = gui.nametowidget('custom_entrybox').get() if (
                gui.nametowidget('custom_entrybox').get() != '') else 'none'
        command.append('custom')
        command.append(str(custom_results_dir))
        command.append(str(custom_import_dir))
        if not queue_contains(command, search_queue):
            modes = (' ' + str(program_modes).replace("'", '')) if len(program_modes) > 0 else ''
            print(f'Enqueueing {tl.short_to_full_date(first_date)}'
                  f'{" - " + tl.short_to_full_date(second_date) if first_date != second_date else ""}'
                  f' on {detector.upper()}{modes}.')
            search_queue.put_nowait(command)


# Starts the search program when the start button is clicked
def start(gui, event, search_queue, program_modes):
    enqueue(gui, search_queue, program_modes)
    if not search_queue.empty():
        # Runs the search manager in a different thread to prevent the GUI from locking up
        search_thread = threading.Thread(target=run, args=(0, gui, event, search_queue))
        search_thread.start()


# Stops the search script when the stop button is clicked
def stop(event):
    if event.running:
        print('Halting program execution...')
        event.set()


# Here to redirect stdout and stderr from the search program
def program_wrapper(write, first_date, second_date, unit, mode_info):
    PipeStdout(write)
    try:
        program(first_date, second_date, unit, mode_info)
    except Exception as ex:
        print('Search program terminated with the following error or warning:\n')
        # Removing the top layer of the traceback (which is just this function) and printing the remainder
        count = len(traceback.extract_tb(ex.__traceback__)) - 1
        print(traceback.format_exc(limit=-count))


# Runs and manages the search program in another process
def run(arg, gui, event, search_queue):  # The useless arg is unfortunately necessary or threading complains
    if not search_queue.empty():
        disable_elements(gui)
        event.running = True

    while not search_queue.empty():
        info = search_queue.get_nowait()
        first_date = info[0]
        second_date = info[1]
        unit = info[2]
        mode_info = info[3:]

        # Prints feedback about what date and modes were selected
        first_date_sep = f'{first_date[0:2]}/{first_date[2:4]}/{first_date[4:]}'
        second_date_sep = f'{second_date[0:2]}/{second_date[2:4]}/{second_date[4:]}'
        print(f'Running search for {first_date_sep}{(" - " + second_date_sep) if first_date != second_date else ""} '
              f'on {unit}.')
        if len(mode_info) > 4:  # one for gui, one for custom, the last two for custom import/export locations
            print(f'This search will be run with the following modes: {", ".join(mode_info[1:-3])}.')

        # Runs the search program in a separate process and manages it
        read, write = multiprocessing.Pipe()
        process = multiprocessing.Process(target=program_wrapper,
                                          args=(write, first_date, second_date, unit, mode_info))
        process.start()
        while process.is_alive() and not event.is_set():
            # Prints the processes' piped stdout (which in turn gets printed to the text box in the gui)
            if read.poll():
                print(read.recv(), end='')

        if process.is_alive():  # This will be executed if the user presses the stop button
            process.terminate()
            break

    event.clear()
    event.running = False
    print('Search Concluded.')
    enable_elements(gui)


# Updates the search queue counter
def update_counter(gui, search_queue):
    gui.nametowidget('enqueue_counter')['text'] = f'Searches\nEnqueued:\n{search_queue.qsize()}'
    milliseconds = 20
    gui.after(milliseconds, update_counter, gui, search_queue)


# Resets all the text entry boxes and tick boxes, as well as the big text box, when the reset button is clicked.
# Also empties the search queue
def reset(gui, event, search_queue, program_modes, variables):
    while not search_queue.empty():  # Emptying the search queue
        search_queue.get_nowait()

    stop(event)

    text_box = gui.nametowidget('text_box')
    text_box['state'] = tk.NORMAL
    text_box.delete('1.0', 'end')
    text_box['state'] = tk.DISABLED

    date_one = gui.nametowidget('date_one')
    date_two = gui.nametowidget('date_two')
    date_one.delete(0, 'end')
    date_one.insert(0, 'yymmdd')
    date_two.delete(0, 'end')
    date_two.insert(0, 'yymmdd')

    gui.nametowidget('detector_entrybox').delete(0, 'end')
    gui.nametowidget('results_entrybox').delete(0, 'end')
    gui.nametowidget('custom_entrybox').delete(0, 'end')

    while program_modes:  # Emptying the modes list
        program_modes.pop()

    for variable in variables:  # Unchecking the checkboxes
        variable.set(0)


def main():
    program_modes = []
    search_queue = Queue()  # Threadsafe queue that holds all the enqueued days
    variables = []  # Checkbox on/off variables
    event = Communicator()

    # General GUI
    gui = tk.Tk()
    gui.title("TGF Search")
    gui.geometry("1080x720")
    gui.resizable(False, False)

    # Making the text box
    text_box = tk.Text(gui, height=30, width=100, name='text_box')
    text_box['state'] = tk.DISABLED
    text_box_label = tk.Label(gui, text='Search Output', anchor=tk.CENTER)

    # Redirecting stdout to the GUI text box
    def redirector(inputStr):
        text_box['state'] = tk.NORMAL
        text_box.insert("end", inputStr)
        text_box.yview(tk.END)
        text_box['state'] = tk.DISABLED

    sys.stdout.write = redirector

    # Packing the text box and label in
    text_box.pack()
    text_box_label.pack()

    # Start and Stop buttons
    start_button = tk.Button(gui, height=3, width=10, text='Start', bg='white', name='start_button',
                             command=lambda: start(gui, event, search_queue, program_modes))
    stop_button = tk.Button(gui, height=3, width=10, text='Stop', bg='white', name='stop_button',
                            command=lambda: stop(event))

    start_button.place(x=430, y=510)
    stop_button.place(x=570, y=510)

    # Input/queue reset button
    reset_button = tk.Button(gui, height=4, width=15, text='Reset/\nClear Queue', bg='white', name='reset_button',
                             command=lambda: reset(gui, event, search_queue, program_modes, variables))

    reset_button.place(x=810, y=510)

    # Clear text box button
    clear_button = tk.Button(gui, height=3, width=8, text='Clear\n Text', bg='white', name='clear_button',
                             command=lambda: clear_box(gui))

    clear_button.place(x=950, y=427)

    # Enqueue button
    enqueue_button = tk.Button(gui, height=3, width=8, text='Enqueue', bg='white', name='enqueue_button',
                               command=lambda: enqueue(gui, search_queue, program_modes))

    enqueue_button.place(x=65, y=427)

    # Enqueue counter
    enqueue_counter = tk.Label(gui, text='', name='enqueue_counter')

    enqueue_counter.place(x=65, y=370)

    # Making and placing date entry boxes
    date_one_label = tk.Label(gui, text='Date One:')
    date_one = tk.Entry(gui, width=15, borderwidth=5, name='date_one')
    date_one.insert(0, 'yymmdd')
    date_one.bind("<FocusIn>", lambda e: clear_ghost_text(date_one, 'yymmdd'))

    date_two_label = tk.Label(gui, text='Date Two: ')
    date_two = tk.Entry(gui, width=15, borderwidth=5, name='date_two')
    date_two.insert(0, 'yymmdd')
    date_two.bind("<FocusIn>", lambda e: clear_ghost_text(date_two, 'yymmdd'))

    date_one_label.place(x=150, y=510)
    date_one.place(x=220, y=510)
    date_two_label.place(x=150, y=550)
    date_two.place(x=220, y=550)

    # Making and placing regular mode checkboxes
    regular_checkbox_label = tk.Label(gui, text='Modes:')

    ascb = tk.IntVar()
    allscints_cb = tk.Checkbutton(gui, text='allscints', variable=ascb, onvalue=1, offvalue=0, name='allscints',
                                  command=lambda: tick_untick(ascb, program_modes, 'allscints'))
    variables.append(ascb)

    tcb = tk.IntVar()
    template_cb = tk.Checkbutton(gui, text='template', variable=tcb, onvalue=1, offvalue=0, name='template',
                                 command=lambda: tick_untick(tcb, program_modes, 'template'))
    variables.append(tcb)

    acb = tk.IntVar()
    aircraft_cb = tk.Checkbutton(gui, text='aircraft', variable=acb, onvalue=1, offvalue=0, name='aircraft',
                                 command=lambda: tick_untick(acb, program_modes, 'aircraft'))
    variables.append(acb)

    combob = tk.IntVar()
    combo_cb = tk.Checkbutton(gui, text='combo', variable=combob, onvalue=1, offvalue=0, name='combo',
                              command=lambda: tick_untick(combob, program_modes, 'combo'))
    variables.append(combob)

    regular_checkbox_label.place(x=150, y=600)
    allscints_cb.place(x=220, y=600)
    template_cb.place(x=300, y=600)
    aircraft_cb.place(x=380, y=600)
    combo_cb.place(x=460, y=600)

    # Making and placing developer mode checkboxes
    dev_checkbox_label = tk.Label(gui, text='Dev Modes:')

    sccb = tk.IntVar()
    skcali_cb = tk.Checkbutton(gui, text='skcali', variable=sccb, onvalue=1, offvalue=0, name='skcali',
                               command=lambda: tick_untick(sccb, program_modes, 'skcali'))
    variables.append(sccb)

    sscb = tk.IntVar()
    skshort_cb = tk.Checkbutton(gui, text='skshort', variable=sscb, onvalue=1, offvalue=0, name='skshort',
                                command=lambda: tick_untick(sscb, program_modes, 'skshort'))
    variables.append(sscb)

    sgcb = tk.IntVar()
    skglow_cb = tk.Checkbutton(gui, text='skglow', variable=sgcb, onvalue=1, offvalue=0, name='skglow',
                               command=lambda: tick_untick(sgcb, program_modes, 'skglow'))
    variables.append(sgcb)

    pcb = tk.IntVar()
    pickle_cb = tk.Checkbutton(gui, text='pickle', variable=pcb, onvalue=1, offvalue=0, name='pickle',
                               command=lambda: tick_untick(pcb, program_modes, 'pickle'))
    variables.append(pcb)

    dev_checkbox_label.place(x=150, y=630)
    skcali_cb.place(x=220, y=630)
    skshort_cb.place(x=300, y=630)
    skglow_cb.place(x=380, y=630)
    pickle_cb.place(x=460, y=630)

    # Making and placing detector entry box
    detector_label = tk.Label(gui, text='Detector:')
    detector_entrybox = tk.Entry(gui, width=10, borderwidth=5, name='detector_entrybox')
    supported_label = tk.Label(gui, text='Supported\n Detectors:\nTHOR(1-6), GODOT,\nSANTIS, CROATIA')

    detector_label.place(x=660, y=615)
    detector_entrybox.place(x=720, y=615)
    supported_label.place(x=800, y=595)

    # Making and placing custom export location entry box and directory dialogue box button
    results_label = tk.Label(gui, text='Export Location:')
    results_entrybox = tk.Entry(gui, width=30, borderwidth=5, name='results_entrybox')
    results_button = tk.Button(gui, width=6, height=2, text='Browse', name='results_button',
                               command=lambda: select_dir(results_entrybox), bg='white')

    results_label.place(x=575, y=680)
    results_entrybox.place(x=675, y=680)
    results_button.place(x=880, y=673)

    # Making and placing custom import location entry box and directory dialogue box button
    custom_label = tk.Label(gui, text='Import Location:')
    custom_entrybox = tk.Entry(gui, width=30, borderwidth=5, name='custom_entrybox')
    custom_button = tk.Button(gui, width=6, height=2, text='Browse', name='custom_buttom',
                              command=lambda: select_dir(custom_entrybox), bg='white')

    custom_label.place(x=150, y=680)
    custom_entrybox.place(x=250, y=680)
    custom_button.place(x=455, y=673)

    # Gui async loop
    update_counter(gui, search_queue)
    tk.mainloop()


if __name__ == '__main__':
    main()
