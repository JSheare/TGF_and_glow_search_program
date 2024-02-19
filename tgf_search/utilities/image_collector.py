"""A script that moves all scatter plots and histograms generated by the search to centralized locations in the pwd."""
import glob
import os
import shutil
import sys

import tgf_search.tools as tl


def main():
    try:
        first_date = str(sys.argv[1])
        second_date = str(sys.argv[2])
        unit = str(sys.argv[3])
    except IndexError:
        print('Please provide a first date, a second date, and a unit name.')
        exit()

    # For the scatter plots
    s_path = f'{os.getcwd()}/Collected Images/Scatter Plots/'
    tl.make_path(s_path)

    # For the histograms
    h_path = f'{os.getcwd()}/Collected Images/Histograms/'
    tl.make_path(h_path)

    detector_path = f"{os.getcwd()}/Results/{unit}/"

    requested_dates = tl.make_date_list(first_date, second_date)
    for date_str in requested_dates:
        path = f'{detector_path}/{date_str}'

        # Scatter plots:
        scatter_plot_list = glob.glob(f'{path}/scatterplots/*event*.png')
        for s_file in scatter_plot_list:
            s_filename = s_file.replace(path, '')[14:]
            shutil.copyfile(s_file, f'{s_path}{s_filename}')

        # Histograms
        histogram_list = glob.glob(f'{path}/*histogram.png')
        for h_file in histogram_list:
            h_filename = h_file.replace(path, '')[1:]
            shutil.copyfile(h_file, f'{h_path}{h_filename}')


if __name__ == '__main__':
    main()
