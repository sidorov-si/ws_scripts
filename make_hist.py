#!/usr/bin/env python
"""
Take a column of numbers as input and plot a histogram 
with a specified bin size.

Usage:
  make_hist.py -i <input_file> -o <output_png_file> -b <bin_size> [--xlabel <histogram_xlabel> --ylabel <histogram_ylabel> --header <histogram_header> --max-value <max_hist_value>]

Options:
  -h --help                     Show this screen.
  --version                     Show version.
  -i <input_filename>           Input file with a column of numbers.
  -o <output_png_filename>      Output file with a histogram.
  -b <bin_size>                 Bin size (interger greater than 0).
  --xlabel <histogram_xlabel>   X label for a histogram. Default: "Count".
  --ylabel <histogram_ylabel>   Y label for a histogram. Default: "Value range".
  --header <histogram_header>   Header for histogram.
  --max-value <max_hist_value>  Threshold for a right border of a range: only ranges with the right border less or equal to the threshold will be plotted.
"""


import sys

print

modules = ["docopt", "os"]
exit_flag = False
for module in modules:
    try:
        __import__(module)
    except ImportError:
        exit_flag = True
        sys.stderr.write("Error: Python module " + module + " is not installed.\n")

if exit_flag:
    sys.stderr.write("You can install these modules with a command: pip install <module>\n")
    sys.stderr.write("(Administrator privileges may be required.)\n")
    sys.exit(1)


from docopt import docopt
from sys import stdout
from os.path import exists
from os.path import isfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pylab import *


def autolabel(rects, max_value):
    # attach some text labels
    for rect in rects:
        width = rect.get_width()
        text_x_coord = 1.01 * width if max_value <= 50 else width + 5
        text(text_x_coord, rect.get_y() + 0.2, '%d'%int(width))


def plot_hist(bin_dict, bin_size, hist_xlabel, hist_ylabel, plot_header, \
              max_hist_value, output_filename):
    sys.stderr.write("Generate histogram...")
    sys.stderr.flush()
    n = len(bin_dict)
    values = []
    is_first_key = True
    previous_key = 0
    last_key = 0
    bin_dict_len = len(bin_dict)
    for key, value in sorted(bin_dict.iteritems()):
        if key > previous_key + 1:
            for i in range(previous_key + 1, key):
                if max_hist_value > 0 and \
                   (i + 1) * bin_size <= max_hist_value or \
                   max_hist_value < 0:
                    values.append(0)
        previous_key = key
        values.append(value)
        last_key = key
    n = len(values)
    pos = [i for i in np.arange(n)]
    if max_hist_value > 0 and max_hist_value <= last_key * bin_size:
        pos = pos[:int(max_hist_value // bin_size)]
        values = values[:int(max_hist_value // bin_size)]
        bar_number = len(values)
    else:
        bar_number = n
    figure(1, figsize=(10, bar_number // 5))
    bars = barh(pos, values, height=1.0, align='center')
    title(plot_header)
    max_value = max(values)
    autolabel(bars, max_value)
    yticks_list = ["(" + str((i - 1) * bin_size) + ", " + \
                   str(i * bin_size) + "]" for i in range(1, bar_number + 1)]
    yticks(pos, yticks_list)
    xlabel(hist_xlabel)
    ylabel(hist_ylabel)
    ylim(ymin = -1, ymax = bar_number)
    xlim(xmin = 0, xmax = max(values) + 0.2 * max(values))
    savefig(output_filename)
    sys.stderr.write("Done.\n\n")
    sys.stderr.flush()


def print_hist(bin_dict, bin_size, max_hist_value):
    previous_key = 0
    for key, value in sorted(bin_dict.items()):
        if int(key) > previous_key + 1:
            for i in range(previous_key + 1, key):
                if max_hist_value > 0 and \
                   i * bin_size <= max_hist_value or \
                   max_hist_value < 0:
                    print "(" + str((int(i) - 1) * bin_size) + ", " + \
                        str(int(i) * bin_size) + "]\t", 0
        if max_hist_value > 0 and \
           key * bin_size <= max_hist_value or \
           max_hist_value < 0:
            print "(" + str((int(key) - 1) * bin_size) + ", " + \
                  str(int(key) * bin_size) + "]\t", value
        previous_key = key
    sys.stderr.write("\n")
    sys.stderr.flush()


def bin_values(input_filename, bin_size, bin_size_type, bin_dict):
    with open(input_filename, 'r') as infile:
        for index, line in enumerate(infile):
            if bin_size_type == "float":
                value = float(line.strip())
            else:
                value = int(line.strip())
            bin_number = int(value // bin_size + 1)
            if bin_number in bin_dict:
                bin_dict[bin_number] += 1
            else:
                bin_dict[bin_number] = 1
            if index % 1000 == 0 and index != 0:
                sys.stderr.write("Processed " + str(index) + " records.\n")
                sys.stderr.flush()
        sys.stderr.write("\n")
        sys.stderr.flush()


if __name__ == '__main__':
    arguments = docopt(__doc__, version='plot_hist 0.6')
    input_filename = arguments["-i"]
    if not exists(input_filename):
        sys.stderr.write("Error: Can't find input file: no such file '" + \
                         input_filename + "'. Exit.\n")
        sys.exit(1)
    if not isfile(input_filename):
        sys.stderr.write("Error: input file must be a regular file. " + \
                         "Something else given. Exit.\n")
        sys.exit(1)
    output_filename = arguments["-o"]
    bin_size_type = "int"
    if "." in arguments["-b"]:
        bin_size = float(arguments["-b"])
        bin_size_type = "float"
    else:
        bin_size = int(arguments["-b"])
    if bin_size < 0:
        sys.stderr.write("Error: bin size must be > 0.\n")
        sys.exit(1)
    if arguments["--max-value"] != None:
        if "." in arguments["--max-value"]:
            max_hist_value = float(arguments["--max-value"])
        else:
            max_hist_value = int(arguments["--max-value"])
        if max_hist_value < bin_size:
            sys.stderr.write("Error: range max value must be >= bin size. Exit.\n")
            sys.exit(1)
    else:
        max_hist_value = -1
    if arguments["--xlabel"] != None:
        hist_xlabel = arguments["--xlabel"]
    else:
        hist_xlabel = "Count"
    if arguments["--ylabel"] != None:
        hist_ylabel = arguments["--ylabel"]
    else:
        hist_ylabel = "Value range"
    if arguments["--header"] != None:
        plot_header = arguments["--header"]
    else:
        plot_header = "Histogram"

    bin_dict = {}

    bin_values(input_filename, bin_size, bin_size_type, bin_dict)
    print_hist(bin_dict, bin_size, max_hist_value)
    plot_hist(bin_dict, bin_size, hist_xlabel, hist_ylabel, plot_header, \
              max_hist_value, output_filename)

