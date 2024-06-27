"""A script that moves all scatter plots and histograms generated by the search to centralized locations in the pwd."""
import glob
import os
import shutil
import sys

# Adds grandparent directory to sys.path. Necessary to make the import below work when running this file as a script
grandparent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
if grandparent_dir not in sys.path:
    sys.path.append(grandparent_dir)

import tgfsearch.tools as tl


def main():
    if len(sys.argv) >= 4:
        first_date = str(sys.argv[1])
        second_date = str(sys.argv[2])
        unit = str(sys.argv[3])
    else:
        print('Please provide a first date, a second date, and a unit name.')
        exit()

    # For the traces
    t_path = f'{os.getcwd()}/Collected Images/Traces'
    tl.make_path(t_path)

    # For the scatter plots
    s_path = f'{os.getcwd()}/Collected Images/Scatter Plots'
    tl.make_path(s_path)

    # For the histograms
    h_path = f'{os.getcwd()}/Collected Images/Histograms'
    tl.make_path(h_path)

    detector_path = f'{os.getcwd()}/Results/{unit}'

    requested_dates = tl.make_date_list(first_date, second_date)
    for date_str in requested_dates:
        path = f'{detector_path}/{date_str}'

        # Traces
        trace_list = glob.glob(f'{path}/traces/*xtr*.png')
        for t_file in trace_list:
            t_filename = t_file.replace(path, '')[8:]
            shutil.copyfile(t_file, f'{t_path}/{t_filename}')

        # Scatter plots:
        scatter_plot_list = glob.glob(f'{path}/scatter_plots/*event*.png')
        for s_file in scatter_plot_list:
            s_filename = s_file.replace(path, '')[14:]
            shutil.copyfile(s_file, f'{s_path}/{s_filename}')

        # Histograms
        histogram_list = glob.glob(f'{path}/*histogram.png')
        for h_file in histogram_list:
            h_filename = h_file.replace(path, '')[1:]
            shutil.copyfile(h_file, f'{h_path}/{h_filename}')


if __name__ == '__main__':
    main()
