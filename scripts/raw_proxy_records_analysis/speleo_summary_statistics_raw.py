# -*- coding: utf-8 -*-
#
# Copyright (C) 2008-2014 Jonathan F. Donges
# Author: Jonathan F. Donges <donges@pik-potsdam.de>
# URL: <http://www.pik-potsdam.de/members/donges/software>

"""
Computes and plots summary statistics from recurrence analysis of several
paleoclimate records.

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

#  Import progress bar for easy progress bar handling
import progressbar

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
RCC_EPISODES = np.array([[8., 9.], [6., 5.], [3.8, 4.2], [2.5, 3.5], [1., 1.2], [0.15, 0.6]]) # See Mayewski et al., Quat. Sci. Rev. (2004)

#  Settings for computation of summary statistics
INTERP_TIME_RES = 10. # years [a]
SIGN_WINDOW_SIZE = 100. # years [a]
N_SURROGATES_SUMM_STAT = 1000 # Number of surrogates for Monte Carlo significance testing of summary statistics


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
#  Functions
#

def summary_statistics(is_significant, time_interp, n_records_at_t, n_window):
    """Compute summary statistics."""
    #  Initialize
    significant_number = np.zeros(len(time_interp))
    significant_relative = np.zeros(len(time_interp))
    max_records_around_t = np.zeros(len(time_interp))

    for j in xrange(len(time_interp)):
        #  Set window and treat boundaries
        start = max(0, j - n_window)
        end = min(j + n_window, len(time_interp))

        #  If any significant values inside a certain window, add an event / a transition
        significant_number[j] = np.any(is_significant[:, start:end], axis=1).sum()

        #  Get max number of records available in window around time j
        max_records_around_t[j] = n_records_at_t[start:end].max()

    #significant_relative[measure] = significant_number[measure] / n_records_at_t
    return (max_records_around_t, significant_number / max_records_around_t)

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
#  Prepare summary statistics of significant transitions/events
#

#  Initialize
significant_number = {}
significant_relative = {}
event_length_dist = {}

#  Get significance window size
n_window = int(SIGN_WINDOW_SIZE / INTERP_TIME_RES)

#  Create regular time axis for interpolation of results
time_interp = np.arange(min_age, max_age, INTERP_TIME_RES)

for measure in symbols.keys():
    #  Initialize
    results_interp = np.zeros((n_time_series, len(time_interp)))
    is_significant = np.zeros((n_time_series, len(time_interp)), dtype=bool)
    event_lengths = np.zeros((n_time_series, len(time_interp)), dtype="int32")

    for i in xrange(n_time_series):
        #  Interpolate results piecewise linearly
        results_interp[i,:] = np.interp(time_interp, time_list[i][step_sequence[i]], results[measure][i], left=np.nan, right=np.nan)
        #  Determine at which time steps the values are outside the confidence band indicating events/transitions
        is_significant[i,:] = np.logical_or(results_interp[i,:] < q_left[measure][i], results_interp[i,:] > q_right[measure][i])

        #  Get frequency distribution of significant event length
        count = 0
        for j in xrange(len(time_interp)):
            if is_significant[i,j]:
                count += 1
                #  Treat boundary
                if j == len(time_interp) - 1:
                    event_lengths[i,count] += 1

            elif 0 < j < len(time_interp) - 1:
                event_lengths[i,count] += 1
                count = 0

    #  Add this to dictionary for each measure
    event_length_dist[measure] = event_lengths

    #  Get number of records available at each time
    n_records_at_t = n_time_series - np.isnan(results_interp).sum(axis=0)

    #  Compute summary statistics
    (max_records_around_t, significant_relative[measure]) = summary_statistics(is_significant, time_interp, n_records_at_t, n_window)

    #  Store results
    np.savetxt("time_axis.txt", time_interp)
    np.savetxt("max_records_around_t.txt", max_records_around_t)
    np.savetxt("summary_stat_" + measure + ".txt", significant_relative[measure])

#
#  Significance testing for summary statistics
#

#  Initialize
q_sum_stat_lower = {}
q_sum_stat_upper = {}

#  Create surrogate event series for each record and measure
for measure in symbols.keys():
    print "Generate surrogates for measure:", measure

    #  Initialize progress bar
    progress = progressbar.ProgressBar().start()

    #  Initialize array for holding generated surrogates
    event_surrogates = np.zeros((n_time_series, N_SURROGATES_SUMM_STAT, len(time_interp)), dtype=bool)

    for i in xrange(n_time_series):
        #  Get indices corresponding to time boundaries of present record
        min_time = time_list[i][step_sequence[i]][0]
        max_time = time_list[i][step_sequence[i]][-1]

        min_time_index = min(range(len(time_interp)), key=lambda i: abs(time_interp[i] - min_time))
        max_time_index = min(range(len(time_interp)), key=lambda i: abs(time_interp[i] - max_time))

        #  Repeat until desired number of valid surrogates is generated
        surr_count = 0
        while surr_count < N_SURROGATES_SUMM_STAT:
            #  Create surrogate
            #surrogate_valid = True

            #  Loop over event lengths
            for j in xrange(1,len(time_interp)):
                #  Loop over all events of a given length j
                for k in xrange(event_length_dist[measure][i,j]):
                    #  Randomly draw index where event of length j begins and make sure that event cannot
                    #  exceed time span of original record
                    index = int(np.random.rand() * (max_time_index - min_time_index - j) + min_time_index)

                    #  Check for overlaps with other events
                    #if event_surrogates[i, surr_count, index:index + j].any():
                    #    surrogate_valid = False
                    #else:
                        #  Write event of length j
                    event_surrogates[i, surr_count, index:index + j] = True

            #  Check surrogate
            #if surrogate_valid:
            surr_count += 1

        #  Update progress bar every step
        progress.update(int(100 * i / float(n_time_series)))

    #  Terminate progress bar
    progress.finish()

    #  Initialize
    summary_stat = np.zeros((N_SURROGATES_SUMM_STAT, len(time_interp)))

    #  Initialize progress bar
    progress = progressbar.ProgressBar().start()

    #  Compute summary statistics for each ensemble member
    for j in xrange(N_SURROGATES_SUMM_STAT):
        summary_stat[j,:] = summary_statistics(event_surrogates[:,j,:], time_interp, n_records_at_t, n_window)[1]

        #  Update progress bar every step
        progress.update(int(100 * j / float(N_SURROGATES_SUMM_STAT)))

    #  Terminate progress bar
    progress.finish()

    #  Get confidence bounds for each time step
    q_elements = int(Q_QUANTILE * N_SURROGATES_SUMM_STAT)

    summary_stat.sort(axis=0)
    q_sum_stat_upper[measure] = summary_stat[-q_elements, :]

    #  Store results of significance test
    np.savetxt("summary_stat_ensemble_" + measure + ".txt", summary_stat)

#
#  Rescale time axis
#

for time in time_list:
    time /= 1000.

min_age /= 1000.
max_age /= 1000.

#
#  Plot summary statistics
#

#  Set plotting parameters (for Clim. Past paper)
params = { 'figure.figsize': (8.268, 8.),
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
          'xtick.minor.size': 3,      # minor tick size in points
          'xtick.major.width': 1.,    # major tick width in points
          'xtick.minor.width': 1.,    # minor tick width in points
          'ytick.minor.size': 3,      # minor tick size in points
          'ytick.major.width': 1.,    # major tick width in points
          'ytick.minor.width': 1.,    # minor tick width in points
          'tick.labelsize': 'small'
           }

#pylab.rcParams.update(params)

#  Create new figure
fig = pylab.figure()

#  Plot maximum number of records available around a certain time
ax = fig.add_subplot(len(symbols) + 1, 1, 1)

pylab.plot(time_interp / 1000., max_records_around_t, "#DB232E")
pylab.fill_between(time_interp / 1000., 0, max_records_around_t, color="#F1B08E")

pylab.ylabel(r"$N_r$", fontsize=16)

#pylab.xlim(min_age, max_age)
pylab.xlim(0, 10)
pylab.ylim(0, 10)

#  Add figure label
ax.annotate(figure_labels[0], xy=(-0.14, 1.1),
        xycoords='axes fraction',
        horizontalalignment='left', verticalalignment='top',
        fontname='Helvetica', fontsize=14, fontweight='bold')

#  Set time axis invisible
#ax.axes.get_xaxis().set_visible(False)

#  Enable minor ticks
ax.minorticks_on()

#  Add time labels and axis label to top subplot
#ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")
ax.tick_params(labelbottom='off',labeltop='on')

pylab.xlabel(AGE_LABEL)

#  Set colors for results subplots
colors_line = ["#3A89CB", "#9CD39A"]
colors_fill = ["#DDEAED", "#DBEDDA"]

#  Plot windowed measures
j = 1
for measure in symbols.keys():
    j += 1
    ax = fig.add_subplot(len(symbols) + 1, 1, j)

    #  Prepare lower bound for fill between to highlight episodes with significant excursions from "normal" climate state
    is_significant = q_sum_stat_upper[measure] < significant_relative[measure]

    #  Plot results
    pylab.plot(time_interp / 1000., significant_relative[measure], colors_line[j - 2])
    #  Fill insignificant episodes
    pylab.fill_between(time_interp / 1000., 0, significant_relative[measure], color=colors_fill[j - 2])
    #  Fill significant episodes
    pylab.fill_between(time_interp / 1000., q_sum_stat_upper[measure], significant_relative[measure], where=is_significant, color=colors_line[j - 2])

    #  Plot confidence band for summary statistics
    #pylab.fill_between(time_interp / 1000., q_sum_stat_lower[measure], q_sum_stat_upper[measure], color="0.75")
    #pylab.plot(time_interp / 1000., q_sum_stat_lower[measure], "r-")
    pylab.plot(time_interp / 1000., q_sum_stat_upper[measure], "-", color="grey")

    pylab.ylabel("$n_" + symbols[measure][1:-1] + "$", fontsize=16)

    #  Plot additional information
    #  Plot RCC episodes
    for k in xrange(RCC_EPISODES.shape[0]):
        pylab.axvspan(RCC_EPISODES[k,0], RCC_EPISODES[k,1], fill=True, edgecolor="0.8", color="0.8", zorder=0)
    #  Plot Bond events
    pylab.scatter(BOND_EVENTS, 0.9 * np.ones(len(BOND_EVENTS)), c='k', marker='o', s=80)

    #pylab.xlim(min_age, max_age)
    pylab.xlim(0, 10)
    pylab.ylim(0, 1)

    #  Add figure label
    ax.annotate(figure_labels[j - 1], xy=(-0.14, 1.1),
                xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontname='Helvetica', fontsize=14, fontweight='bold')

    #  Enable minor ticks
    ax.minorticks_on()

    #  Render time tick labels invisible
    ax.set_xticklabels([])

    #  Set number of ticks and tick labels
    ax.locator_params(axis="y", tight=True, nbins=5)


#  Enable minor ticks
ax.minorticks_on()

pylab.xlabel(AGE_LABEL)

#  Save figure
pylab.savefig("summary_statistics.pdf")

#  Show all figures
pylab.show()
