import tkinter as tk
from tkinter import filedialog
import sys as sys
import subprocess as subprocess
import os as os
import signal as signal
import threading as threading


# Need this to make sure that the search queue isn't emptied when the user presses the stop button
class lockableQueue:
    def __init__(self):
        self.queue = []
        self.locked = False

    def __bool__(self):
        return len(self.queue) > 0

    def __contains__(self, item):
        return item in self.queue

    def __len__(self):
        return len(self.queue)

    def push(self, item):
        if not self.locked:
            self.queue.insert(0, item)

    def pop(self):
        if not self.locked:
            return self.queue.pop()

    def lock(self):
        self.locked = True

    def unlock(self):
        self.locked = False


# Creates a directory selection dialogue box and then puts the selected directory in the specified text entry box
def select_dir(entry_box):
    directory = filedialog.askdirectory(initialdir='/')
    entry_box.delete(0, 'end')
    entry_box.insert(0, directory)


# Clears the sample text from the date entry boxes when they are clicked
def ghost_text_clear(entrybox, ghost_text):
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
def clear_box():
    text_box['state'] = tk.NORMAL
    text_box.delete('1.0', 'end')
    text_box['state'] = tk.DISABLED


# Disables all checkboxes/buttons
def disable_elements():
    global regular_checkbox_list, dev_checkbox_list
    start_button['state'] = tk.DISABLED
    enqueue_button['state'] = tk.DISABLED
    date_one['state'] = tk.DISABLED
    date_two['state'] = tk.DISABLED
    detector_entrybox['state'] = tk.DISABLED
    custom_entrybox['state'] = tk.DISABLED
    results_entrybox['state'] = tk.DISABLED
    for box in regular_checkbox_list + dev_checkbox_list:
        box['state'] = tk.DISABLED


# Enables all checkboxes/buttons
def enable_elements():
    global regular_checkbox_list, dev_checkbox_list
    start_button['state'] = tk.NORMAL
    enqueue_button['state'] = tk.NORMAL
    date_one['state'] = tk.NORMAL
    date_two['state'] = tk.NORMAL
    detector_entrybox['state'] = tk.NORMAL
    custom_entrybox['state'] = tk.NORMAL
    results_entrybox['state'] = tk.NORMAL
    for box in regular_checkbox_list + dev_checkbox_list:
        box['state'] = tk.NORMAL


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


def enqueue():
    global search_queue, program_modes
    first_date = date_one.get()
    second_date = date_two.get()
    detector = detector_entrybox.get()
    # If the search command is valid, constructs the command and adds it to the queue
    if is_valid_search(first_date, second_date, detector):
        command = f'python3 -u search.py {first_date} {second_date} {detector.upper()}'.split()
        for mode in program_modes:
            command.append(mode)

        custom_results_dir = results_entrybox.get() if results_entrybox.get() != '' else 'none'
        custom_import_dir = custom_entrybox.get() if custom_entrybox.get() != '' else 'none'
        command.append('GUI')
        command.append(str(custom_results_dir))
        command.append(str(custom_import_dir))
        if command not in search_queue and not search_queue.locked:
            first_date_sep = f'{first_date[0:2]}/{first_date[2:4]}/{first_date[4:]}'
            second_date_sep = f'{second_date[0:2]}/{second_date[2:4]}/{second_date[4:]}'
            modes = (' ' + str(program_modes).replace("'", '')) if len(program_modes) > 0 else ''
            print(f'Enqueueing {first_date_sep}'
                  f'{" - " + second_date_sep if first_date != second_date else ""}'
                  f' on {detector.upper()}{modes}.')
            search_queue.push(command)


# Starts the search script when the start button is clicked
def start():
    global search_queue, program_modes
    enqueue()
    if search_queue:
        disable_elements()

        # Runs the search script in a different thread to prevent the GUI from locking up
        search_thread = threading.Thread(target=run, args=(0, ))
        search_thread.start()


# Stops the search script from running and unlocks the start button/tick boxes when the stop button is clicked
def stop():
    global pid, search_queue
    # Kills the program
    if pid is not None:
        print('Halting program execution (this might take a second)')
        search_queue.lock()
        os.kill(pid, signal.SIGTERM)

    pid = None
    enable_elements()


# Runs the search script and pipes stdout into the stdout queue
def run(arg):  # The useless arg is unfortunately necessary or threading will complain
    global pid, stdout_queue, search_queue
    while search_queue and not search_queue.locked:
        command = search_queue.pop()
        # Prints feedback about what date and modes were selected
        first_date = command[3]
        second_date = command[4]
        first_date_sep = f'{first_date[0:2]}/{first_date[2:4]}/{first_date[4:]}'
        second_date_sep = f'{second_date[0:2]}/{second_date[2:4]}/{second_date[4:]}'
        stdout_queue.append(f'Running search for {first_date_sep}'
                            f'{(" - " + second_date_sep) if first_date != second_date else ""}'
                            f' on {command[5]}.')

        if len(command) > 9:
            stdout_queue.append(f'This search will be run with the following modes: {", ".join(command[6:-3])}.')

        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pid = process.pid

        while True:
            output = process.stdout.readline()
            if process.poll() is not None:
                break
            if output:
                stdout_queue.append(output.strip().decode('utf-8'))

        pid = None

        out, err = process.communicate()
        if err:
            stdout_queue.append('\n')
            stdout_queue.append('Search script terminated with the following error or warning:')
            stdout_queue.append(err.strip().decode('utf-8'))
            stdout_queue.append('\n')

    stdout_queue.append('Search Concluded.')


# Checks to see if there are items in the queue and prints them, updates the search queue counter, and unlocks the start
# button/tick boxes once the search script has finished executing
def checker():
    global gui, stdout_queue, search_queue
    while len(stdout_queue) > 0:
        text = stdout_queue.pop(0)
        print(text)
        if text == 'Search Concluded.':
            enable_elements()
            search_queue.unlock()

    enqueue_counter['text'] = f'Searches\nEnqueued:\n{len(search_queue)}'

    milliseconds = 20
    gui.after(milliseconds, checker)


# Resets all the text entry boxes and tick boxes, as well as the big text box, when the reset button is clicked.
# Also clears the search queue
def reset():
    global search_queue, program_modes, regular_checkbox_variables, dev_checkbox_variables
    search_queue = lockableQueue()
    stop()

    text_box['state'] = tk.NORMAL
    text_box.delete('1.0', 'end')
    text_box['state'] = tk.DISABLED

    date_one.delete(0, 'end')
    date_one.insert(0, 'yymmdd')
    date_two.delete(0, 'end')
    date_two.insert(0, 'yymmdd')

    detector_entrybox.delete(0, 'end')
    results_entrybox.delete(0, 'end')
    custom_entrybox.delete(0, 'end')

    program_modes = []

    for box in regular_checkbox_variables + dev_checkbox_variables:
        box.set(0)


if __name__ == '__main__':
    program_modes = []
    search_queue = lockableQueue()  # Queue that holds all the enqueued days
    stdout_queue = []  # Queue that holds all the stdout output strings from the search script until they can be printed
    pid = None  # Program identification number for the search script

    # General GUI
    gui = tk.Tk()
    gui.title("TGF Search")
    gui.geometry("1080x720")
    gui.resizable(False, False)

    # Making the text box
    text_box = tk.Text(gui, height=30, width=100)
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
    start_button = tk.Button(gui, height=3, width=10, text='Start',
                             command=lambda: start(), bg='white')
    stop_button = tk.Button(gui, height=3, width=10, text='Stop',
                            command=lambda: stop(), bg='white')

    start_button.place(x=430, y=510)
    stop_button.place(x=570, y=510)

    # Input/queue reset button
    reset_button = tk.Button(gui, height=4, width=15, text='Reset/\nClear Queue',
                             command=lambda: reset(), bg='white')

    reset_button.place(x=810, y=510)

    # Clear text box button
    clear_button = tk.Button(gui, height=3, width=8, text='Clear\n Text',
                             command=lambda: clear_box(), bg='white')

    clear_button.place(x=950, y=427)

    # Enqueue button
    enqueue_button = tk.Button(gui, height=3, width=8, text='Enqueue',
                               command=lambda: enqueue(), bg='white')

    enqueue_button.place(x=65, y=427)

    # Enqueue counter
    enqueue_counter = tk.Label(gui, text='')

    enqueue_counter.place(x=65, y=370)

    # Making and placing date entry boxes
    date_one_label = tk.Label(gui, text='Date One:')
    date_one = tk.Entry(gui, width=15, borderwidth=5)
    date_one.insert(0, 'yymmdd')
    date_one.bind("<FocusIn>", lambda e: ghost_text_clear(date_one, 'yymmdd'))

    date_two_label = tk.Label(gui, text='Date Two: ')
    date_two = tk.Entry(gui, width=15, borderwidth=5)
    date_two.insert(0, 'yymmdd')
    date_two.bind("<FocusIn>", lambda e: ghost_text_clear(date_two, 'yymmdd'))

    date_one_label.place(x=150, y=510)
    date_one.place(x=220, y=510)
    date_two_label.place(x=150, y=550)
    date_two.place(x=220, y=550)

    # Making and placing regular mode checkboxes
    regular_checkbox_label = tk.Label(gui, text='Modes:')

    ascb = tk.IntVar()
    allscints_cb = tk.Checkbutton(gui, text='allscints', variable=ascb, onvalue=1, offvalue=0,
                                  command=lambda: tick_untick(ascb, program_modes, 'allscints'))

    tcb = tk.IntVar()
    template_cb = tk.Checkbutton(gui, text='template', variable=tcb, onvalue=1, offvalue=0,
                                 command=lambda: tick_untick(tcb, program_modes, 'template'))

    acb = tk.IntVar()
    aircraft_cb = tk.Checkbutton(gui, text='aircraft', variable=acb, onvalue=1, offvalue=0,
                                 command=lambda: tick_untick(acb, program_modes, 'aircraft'))

    combob = tk.IntVar()
    combo_cb = tk.Checkbutton(gui, text='combo', variable=combob, onvalue=1, offvalue=0,
                              command=lambda: tick_untick(combob, program_modes, 'combo'))

    regular_checkbox_label.place(x=150, y=600)
    allscints_cb.place(x=220, y=600)
    template_cb.place(x=300, y=600)
    aircraft_cb.place(x=380, y=600)
    combo_cb.place(x=460, y=600)

    regular_checkbox_list = [allscints_cb, template_cb, aircraft_cb, combo_cb]
    regular_checkbox_variables = [ascb, tcb, acb, combob]

    # Making and placing developer mode checkboxes
    dev_checkbox_label = tk.Label(gui, text='Dev Modes:')

    sccb = tk.IntVar()
    skcali_cb = tk.Checkbutton(gui, text='skcali', variable=sccb, onvalue=1, offvalue=0,
                               command=lambda: tick_untick(sccb, program_modes, 'skcali'))

    sscb = tk.IntVar()
    skshort_cb = tk.Checkbutton(gui, text='skshort', variable=sscb, onvalue=1, offvalue=0,
                                command=lambda: tick_untick(sscb, program_modes, 'skshort'))

    sgcb = tk.IntVar()
    skglow_cb = tk.Checkbutton(gui, text='skglow', variable=sgcb, onvalue=1, offvalue=0,
                               command=lambda: tick_untick(sgcb, program_modes, 'skglow'))

    pcb = tk.IntVar()
    pickle_cb = tk.Checkbutton(gui, text='pickle', variable=pcb, onvalue=1, offvalue=0,
                               command=lambda: tick_untick(pcb, program_modes, 'pickle'))

    dev_checkbox_label.place(x=150, y=630)
    skcali_cb.place(x=220, y=630)
    skshort_cb.place(x=300, y=630)
    skglow_cb.place(x=380, y=630)
    pickle_cb.place(x=460, y=630)

    dev_checkbox_list = [skcali_cb, skshort_cb, skglow_cb, pickle_cb]
    dev_checkbox_variables = [sccb, sscb, sgcb, pcb]

    # Making and placing detector entry box
    detector_label = tk.Label(gui, text='Detector:')
    detector_entrybox = tk.Entry(gui, width=10, borderwidth=5)
    supported_label = tk.Label(gui, text='Supported\n Detectors:\nTHOR(1-6), GODOT,\nSANTIS, CROATIA')

    detector_label.place(x=660, y=615)
    detector_entrybox.place(x=720, y=615)
    supported_label.place(x=800, y=595)

    # Making and placing custom export location entry box and directory dialogue box button
    results_label = tk.Label(gui, text='Export Location:')
    results_entrybox = tk.Entry(gui, width=30, borderwidth=5)
    results_button = tk.Button(gui, width=6, height=2, text='Browse',
                               command=lambda: select_dir(results_entrybox), bg='white')

    results_label.place(x=575, y=680)
    results_entrybox.place(x=675, y=680)
    results_button.place(x=880, y=673)

    # Making and placing custom import location entry box and directory dialogue box button
    custom_label = tk.Label(gui, text='Import Location:')
    custom_entrybox = tk.Entry(gui, width=30, borderwidth=5)
    custom_button = tk.Button(gui, width=6, height=2, text='Browse',
                              command=lambda: select_dir(custom_entrybox), bg='white')

    custom_label.place(x=150, y=680)
    custom_entrybox.place(x=250, y=680)
    custom_button.place(x=455, y=673)

    # Gui async loop
    checker()
    tk.mainloop()
