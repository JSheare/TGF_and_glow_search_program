import tkinter as tk
import sys as sys
import subprocess as subprocess
import os as os
import signal as signal
import threading as threading


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


# Stops the search script from running and unlocks the start button/tick boxes when the stop button is clicked
def stop(regular_boxes, dev_boxes):
    global pid
    # Kills the program
    if pid is not None:
        print('Halting program execution (this might take a second)')
        os.kill(pid, signal.SIGTERM)

    pid = None
    # Re-enabling buttons
    start_button['state'] = tk.NORMAL
    for box in regular_boxes + dev_boxes:
        box['state'] = tk.NORMAL


# Checks to see if there are items in the queue and prints them
# Also unlocks the start button/tick boxes once the search script has finished executing
def checker(regular_boxes, dev_boxes):
    global gui, queue
    queue = queue[::-1]
    while len(queue) > 0:
        text = queue.pop()
        print(text)
        if text == 'Search Concluded.':
            # Re-enabling buttons
            start_button['state'] = tk.NORMAL
            for box in regular_boxes + dev_boxes:
                box['state'] = tk.NORMAL

    milliseconds = 20
    gui.after(milliseconds, checker, regular_boxes, dev_boxes)


# Runs the search script and pipes stdout into the queue
def run(command, arg):  # The useless arg is unfortunately necessary or threading will complain
    global pid, queue
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pid = process.pid

    while True:
        output = process.stdout.readline()
        if process.poll() is not None:
            break
        if output:
            queue.append(output.strip().decode('utf-8'))

    pid = None

    out, err = process.communicate()
    if err:
        queue.append('\n')
        queue.append('Search script terminated with the following error:')
        queue.append(err.strip().decode('utf-8'))
        queue.append('\n')

    queue.append('Search Concluded.')


# Resets all the text entry boxes and tick boxes, as well as the big text box, when the reset button is clicked
def reset(modes, reg_variables, dev_variables, regular_boxes, dev_boxes):
    stop(regular_boxes, dev_boxes)

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

    for item in range(len(modes)):
        modes.pop()

    for box in reg_variables + dev_variables:
        box.set(0)


# Starts the search script when the start button is clicked
def start(modes, regular_boxes, dev_boxes):
    first_date = date_one.get()
    second_date = date_two.get()

    # Checks that both dates are digits in the proper format
    if not first_date.isdigit() or not second_date.isdigit() \
            or len(first_date) != 6 or len(second_date) != 6:
        print('Error: not a valid date. Please enter BOTH dates in yymmdd format.')
        return

    # Checks that both dates are sequential
    if int(first_date) > int(second_date):
        print('Error: second date must be AFTER first date.')
        return

    # Checks that a detector has been entered
    detector = detector_entrybox.get()
    if detector == '':
        print('Error: Please enter a detector.')
        return

    # Disabling all checkboxes/buttons
    start_button['state'] = tk.DISABLED
    for box in regular_boxes + dev_boxes:
        box['state'] = tk.DISABLED

    # Assembles the command and running it with search.py
    command = f'python3 -u search.py {first_date} {second_date} {detector.upper()}'

    for mode in modes:
        command += ' ' + mode

    custom_results_dir = results_entrybox.get() if results_entrybox.get() != '' else 'none'
    custom_import_dir = custom_entrybox.get() if custom_entrybox.get() != '' else 'none'
    command += f' GUI {custom_results_dir} {custom_import_dir}'

    # Prints feedback about what date and modes were selected
    first_date_sep = f'{first_date[0:2]}/{first_date[2:4]}/{first_date[4:]}'
    second_date_sep = f'{second_date[0:2]}/{second_date[2:4]}/{second_date[4:]}'
    print(f'Running search for '
          f'{first_date_sep + " - " + second_date_sep if first_date != second_date else first_date_sep}.')

    if len(modes) != 0:
        print(f'This search will be run with the following modes: {", ".join(modes)}.')

    # Runs the search script in a different thread to prevent the GUI from locking up
    search_thread = threading.Thread(target=run, args=(command, 0))
    search_thread.start()


program_modes = []
pid = None  # Program identification number for the search script
queue = []  # Queue that holds all the stdout output strings from the search script until they can be printed

# General GUI
gui = tk.Tk()
gui.title("TGF Search")
gui.geometry("1080x720")
gui.resizable(False, False)

# Making the text box
text_box = tk.Text(gui, height=30, width=100)
text_box['state'] = tk.DISABLED
text_box_label = tk.Label(gui, text="Search Output")


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
                         command=lambda: start(program_modes, regular_checkbox_list, dev_checkbox_list),
                         bg='white')
stop_button = tk.Button(gui, height=3, width=10, text='Stop',
                        command=lambda: stop(regular_checkbox_list, dev_checkbox_list), bg='white')

start_button.place(x=430, y=510)
stop_button.place(x=570, y=510)

# Input reset button
reset_button = tk.Button(gui, height=4, width=15, text='Reset',
                         command=lambda: reset(program_modes, regular_checkbox_variables, dev_checkbox_variables,
                                               regular_checkbox_list, dev_checkbox_list), bg='white')

reset_button.place(x=810, y=510)

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

regular_checkbox_label.place(x=150, y=600)
allscints_cb.place(x=220, y=600)
template_cb.place(x=300, y=600)
aircraft_cb.place(x=380, y=600)

regular_checkbox_list = [allscints_cb, template_cb, aircraft_cb]
regular_checkbox_variables = [ascb, tcb, acb]

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

# Making and placing custom export location entry box
results_label = tk.Label(gui, text='Export Location:')
results_entrybox = tk.Entry(gui, width=40, borderwidth=5)

results_label.place(x=150, y=680)
results_entrybox.place(x=250, y=680)

# Making and placing custom import location entry box
custom_label = tk.Label(gui, text='Import Location:')
custom_entrybox = tk.Entry(gui, width=40, borderwidth=5)

custom_label.place(x=575, y=680)
custom_entrybox.place(x=675, y=680)

# Gui async loop
checker(regular_checkbox_list, dev_checkbox_list)
tk.mainloop()
