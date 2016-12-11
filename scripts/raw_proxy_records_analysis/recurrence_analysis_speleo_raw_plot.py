# -*- coding: utf-8 -*-
#
# Copyright (C) 2008-2014 Jonathan F. Donges
# Author: Jonathan F. Donges <donges@pik-potsdam.de>
# URL: <http://www.pik-potsdam.de/members/donges/software>

"""
Plots results of recurrence analysis of paleoclimate proxy records.

This script provides analyses for this publication:

J.F. Donges, R.V. Donner, N. Marwan, S.F.M. Breitenbach, K. Rehfeld, and J. Kurths,
Nonlinear regime shifts in Holocene Asian monsoon variability: Potential impacts on cultural change and migratory patterns,
Climate of the Past 11, 709-741 (2015),
DOI: 10.5194/cp-11-709-2015
"""

#
#  Imports
#

#  Import sys for command line arguments
import sys

#  Import cPickle for loading and saving data
import cPickle

#  Import np for fast numerics
import numpy as np

#  Import pylab for plotting
import pylab

#
#  Settings
#

#  Settings for statistics

#  q-quantile
Q_QUANTILE = 0.10

#  Settings for plotting

#  List of letters from "a" to "z"
#figure_labels = map(chr, range(97, 123))

#  List of letters from "A" to "Z"
figure_labels = map(chr, range(65, 91))

#  Toggle plotting of confidence band
SHOW_CONFIDENCE_BAND = True

#  Set labels
AGE_LABEL = "Age (ka B.P.)"
DATA_LABEL = r"$\Delta \delta^{18}$O" #  (‰ VPDB)

#  Settings for plotting additional information

#  10.3 ka BP, 9.4 ka BP, 8.2 ka BP, 5.9 ka BP, 4.3–4.4 ka BP, 2.8 ka BP, 1.3–1.1 ka BP
BOND_EVENTS = np.array([10.3, 9.4, 8.2, 5.9, 4.35, 2.8, 1.2, 0.5])
RCC_EPISODES = np.array([[8.,9.],[6.,5.],[3.8,4.2],[2.5,3.5],[1.,1.2],[0.15,0.6]]) # See Mayewski et al., Quat. Sci. Rev. (2004)

#
#  Load results for plotting
#

#  Load from file
filename = sys.argv[1]
print "Loading results from", filename

file = open(filename, 'r')
storage = cPickle.load(file)
file.close()

#  Store parameters
FILENAMES = storage["FILENAMES"]
NAMES = storage["NAMES"]

T_TIME = storage["T_TIME"]
DELTA_TIME = storage["DELTA_TIME"]
DETRENDING_WINDOW_SIZE = storage["DETRENDING_WINDOW_SIZE"]

DIM = storage["DIM"]
TAU = storage["TAU"]
METRIC = storage["METRIC"]
RR = storage["RR"]

N_ENSEMBLE = storage["N_ENSEMBLE"]
SHUFFLE_EMBEDDED = storage["SHUFFLE_EMBEDDED"]

#  Store symbols
symbols = storage["symbols"]

#  Store raw input data
time_list = storage["time_list"]
data_list = storage["data_list"]

#  Store axes
step_sequence = storage["step_sequence"]

#  Store results
results = storage["results"]
surrogate_results = storage["surrogate_results"]

#  Derived properties

#  Get number of time series
n_time_series = len(FILENAMES)

print "Number of realizations for significance test:", N_ENSEMBLE

#
#  Prepare time axis
#

min_age_list = []
max_age_list = []

for i in xrange(n_time_series):
    #  Get minimum and maximum age in time series
    min_age_list.append(time_list[i].min())
    max_age_list.append(time_list[i].max())

#  Get overall minimum and maximum age
min_age = np.array(min_age_list).min()
max_age = np.array(max_age_list).max()

#
#  Calculate means and q-quantiles
#

q_left = {}
q_right = {}
surrogate_mean = {}

for measure in symbols.keys():
    q_left[measure] = []
    q_right[measure] = []
    surrogate_mean[measure] = []

#  Number of elements right and left of the quantile
#  Divide by two, such that the total probability of false rejection of null hypothesis
#  is equal to Q_QUANTILE (two sided test).
q_elements = int(Q_QUANTILE / 2 * N_ENSEMBLE)

for i in xrange(n_time_series):
    for measure in symbols.keys():
        surrogate_results[measure][i].sort()
        #  Get mean
        surrogate_mean[measure].append(surrogate_results[measure][i].mean())
        #  Get q-quantiles
        q_left[measure].append(surrogate_results[measure][i][q_elements - 1])
        q_right[measure].append(surrogate_results[measure][i][-q_elements])

#
#  Print results of significance test
#

for i in xrange(n_time_series):
    for measure in symbols.keys():
        print "\n" + measure + ":"
        print "Mean:", surrogate_mean[measure][i]
        print "Left", Q_QUANTILE / 2, "quantile:", q_left[measure][i]
        print "Right", 1 - Q_QUANTILE / 2, "quantile:", q_right[measure][i]

#
#  Print deviations of window sizes
#

for i in xrange(n_time_series):
    print "\n", NAMES[i], ":"

    average_sampling_time = np.diff(time_list[i]).mean()
    window_size_data_points = int(T_TIME[i] / average_sampling_time)
    step_size_data_points = int(DELTA_TIME / average_sampling_time)

    window_size = np.empty(len(step_sequence[i]))

    for j in xrange(len(step_sequence[i])):
        window_size[j] = (time_list[i][step_size_data_points * j + window_size_data_points]
                          - time_list[i][step_size_data_points * j])

    print "Max. window size:", window_size.max()
    print "Min. window size:", window_size.min()
    print "Mean window size:", window_size.mean()
    print "Max. deviation from prescribed window size [time]:", np.abs(window_size - T_TIME[i]).max()
    print "Max. deviation from prescribed window size [percent]:", 100 * np.abs(window_size - T_TIME[i]).max() / T_TIME[i]

    step_size = np.diff(time_list[i][step_sequence[i]])
    print "Max. step size:", step_size.max()
    print "Min. step size:", step_size.min()
    print "Mean step size:", step_size.mean()
    print "Max. deviation from prescribed step size [time]:", np.abs(step_size - DELTA_TIME).max()
    print "Max. deviation from prescribed step size [percent]:", 100 * np.abs(step_size - DELTA_TIME).max() / DELTA_TIME

#
#  Rescale time axis
#

for time in time_list:
    time /= 1000.

min_age /= 1000.
max_age /= 1000.

#
#  Plotting
#

#  Set plotting parameters
params = { 'figure.figsize': (8.268,11.693),
          'axes.labelsize': 12,
          'text.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'legend.fontsize': 10,
          'title.fontsize': 12,
          'text.usetex': False,
          'font': 'Helvetica',
          'mathtext.bf': 'helvetica:bold',
           'xtick.major.pad': 6,
           'ytick.major.pad': 6,
           'xtick.major.size': 5,
           'ytick.major.size': 5,
           'tick.labelsize': 'small'
           }

#pylab.rcParams.update(params)

#
#  Plot overview figure for each time series
#

for i in xrange(n_time_series):
    #  Create new figure
    fig = pylab.figure()

    #  Set title
    pylab.suptitle(NAMES[i])

    #  Plot original data
    ax = fig.add_subplot(len(symbols) + 1,1,1)

    pylab.plot(time_list[i], data_list[i], "k")
    pylab.ylabel(DATA_LABEL)
    pylab.xlim(time_list[i].min(), time_list[i].max())

    #  Add figure label
    ax.annotate(figure_labels[0], xy=(-0.14, 1.1),
                xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontname='Helvetica', fontsize=14, fontweight='bold')

    #  Plot windowed measures
    j = 1
    for measure in symbols.keys():
        j += 1
        ax = fig.add_subplot(len(symbols) + 1,1,j)
        pylab.plot(time_list[i][step_sequence[i]], results[measure][i], "k")
        pylab.ylabel(symbols[measure])

        if SHOW_CONFIDENCE_BAND:
            #pylab.axhline(y=surrogate_mean[measure][i], linewidth=1, color="r")
            pylab.fill_between(time_list[i][step_sequence[i]], q_left[measure][i], q_right[measure][i], color="0.75")

        pylab.xlim(time_list[i].min(), time_list[i].max())

        #  Add figure label
        ax.annotate(figure_labels[j - 1], xy=(-0.14, 1.1),
                xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontname='Helvetica', fontsize=14, fontweight='bold')

    pylab.xlabel(AGE_LABEL)

    #  Save figure
    pylab.savefig(NAMES[i] + "_overview.pdf")

#
#  Plot comparison figures for single time series / windowed measures
#

#  Comparison of time series
fig = pylab.figure()

for i in xrange(n_time_series):
    ax = fig.add_subplot(n_time_series, 1, i+1)

    pylab.plot(time_list[i], data_list[i], "k")

    pylab.xlim(min_age, max_age)

    #  Plot additional information
    #  Plot RCC episodes
    for k in xrange(RCC_EPISODES.shape[0]):
        pylab.axvspan(RCC_EPISODES[k,0], RCC_EPISODES[k,1], hatch="/", fill=False, edgecolor="0.5")
    #  Plot Bond events
    pylab.scatter(BOND_EVENTS, 0.9 * np.ones(len(BOND_EVENTS)) * data_list[i].max(), c='k', marker='*', s=80)

    #  Add record name
    ax.text(0.99, 0.9, NAMES[i], horizontalalignment='right', verticalalignment='top',
                transform = ax.transAxes)

    #  Add figure label
    ax.annotate(figure_labels[i], xy=(-0.14, 1.1),
                xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontname='Helvetica', fontsize=14, fontweight='bold')

    #  Remove spines
    for loc, spine in ax.spines.items():
        #  For first subplot
        if i == 0:
            if loc in ['left','top']:
                spine.set_position(('outward',10)) # outward by 10 points
            elif loc in ['right','bottom']:
                spine.set_color('none') # don't draw spine
            else:
                raise ValueError('unknown spine location: %s'%loc)

            # turn off ticks where there is no spine
            ax.xaxis.set_ticks_position('top')
            ax.yaxis.set_ticks_position('left')
        #  For every even numbered subplot
        elif i % 2 == 0:
            if loc in ['left','bottom']:
                spine.set_position(('outward',10)) # outward by 10 points
            elif loc in ['right','top']:
                spine.set_color('none') # don't draw spine
            else:
                raise ValueError('unknown spine location: %s'%loc)

            # turn off ticks where there is no spine
            #ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        #  For every uneven numbered subplot
        elif i % 2 == 1:
            if loc in ['right','bottom']:
                spine.set_position(('outward',10)) # outward by 10 points
            elif loc in ['left','top']:
                spine.set_color('none') # don't draw spine
            else:
                raise ValueError('unknown spine location: %s'%loc)

            # turn off ticks where there is no spine
            #ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('right')

    #  Add time ticks and label to uppermost subplot
    if i == 0:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        pylab.xlabel(AGE_LABEL)

    #  Move y-axis tick labels and axis label to right side for every second subplot
    if i % 2 == 1:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

        #for label in ax.yaxis.get_ticklabels():
        #    label.set_horizontalalignment('right')

    #  Render time tick labels invisible
    if 0 < i < n_time_series - 1:
        ax.set_xticklabels([])

    #  Set number of ticks and tick labels
    ax.locator_params(axis="x", tight=True, nbins=9)
    ax.locator_params(axis="y", tight=True, nbins=4)

    pylab.ylabel(DATA_LABEL)

pylab.xlabel(AGE_LABEL)

#  Save figure
pylab.savefig("time_series.pdf")

#  Comparison of windowed measures
for measure in symbols.keys():
    fig = pylab.figure()

    for i in xrange(n_time_series):
        ax = fig.add_subplot(n_time_series, 1, i + 1)

        pylab.plot(time_list[i][step_sequence[i]], results[measure][i], "k")

        pylab.xlim(min_age, max_age)

        if SHOW_CONFIDENCE_BAND:
            #pylab.axhline(y=surrogate_mean[measure][i], linewidth=1, color="r")
            pylab.fill_between(time_list[i][step_sequence[i]], q_left[measure][i], q_right[measure][i], color="0.75")

        #  Plot additional information
        #  Plot RCC episodes
        for k in xrange(RCC_EPISODES.shape[0]):
            pylab.axvspan(RCC_EPISODES[k,0], RCC_EPISODES[k,1], hatch="/", fill=False, edgecolor="0.5")
        #  Plot Bond events
        pylab.scatter(BOND_EVENTS, 0.9 * np.ones(len(BOND_EVENTS)) * results[measure][i].max(), c='k', marker='*', s=80)

        #  Add figure label
        ax.annotate(figure_labels[i], xy=(-0.14, 1.1),
                xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontname='Helvetica', fontsize=14, fontweight='bold')

        #  Add record name
        ax.text(0.99, 0.9, NAMES[i], horizontalalignment='right', verticalalignment='top',
                    transform = ax.transAxes)

        #  Remove spines
        for loc, spine in ax.spines.items():
            #  For first subplot
            if i == 0:
                if loc in ['left','top']:
                    spine.set_position(('outward',10)) # outward by 10 points
                elif loc in ['right','bottom']:
                    spine.set_color('none') # don't draw spine
                else:
                    raise ValueError('unknown spine location: %s'%loc)

                # turn off ticks where there is no spine
                ax.xaxis.set_ticks_position('top')
                ax.yaxis.set_ticks_position('left')
            #  For every even numbered subplot
            elif i % 2 == 0:
                if loc in ['left','bottom']:
                    spine.set_position(('outward',10)) # outward by 10 points
                elif loc in ['right','top']:
                    spine.set_color('none') # don't draw spine
                else:
                    raise ValueError('unknown spine location: %s'%loc)

                # turn off ticks where there is no spine
                #ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
            #  For every uneven numbered subplot
            elif i % 2 == 1:
                if loc in ['right','bottom']:
                    spine.set_position(('outward',10)) # outward by 10 points
                elif loc in ['left','top']:
                    spine.set_color('none') # don't draw spine
                else:
                    raise ValueError('unknown spine location: %s'%loc)

                # turn off ticks where there is no spine
                #ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('right')

        #  Add time ticks and label to uppermost subplot
        if i == 0:
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position("top")
            pylab.xlabel(AGE_LABEL)

        #  Move y-axis tick labels and axis label to right side for every second subplot
        if i % 2 == 1:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")

            #for label in ax.yaxis.get_ticklabels():
            #    label.set_horizontalalignment('right')

        #  Render time tick labels invisible
        if 0 < i < n_time_series - 1:
            ax.set_xticklabels([])

        #  Set number of ticks and tick labels
        ax.locator_params(axis="x", tight=True, nbins=9)
        ax.locator_params(axis="y", tight=True, nbins=4)

        pylab.ylabel(symbols[measure])

    pylab.xlabel(AGE_LABEL)

    #  Save figure
    pylab.savefig("RN_" + measure + ".pdf")

#  Show all figures
#pylab.show()
