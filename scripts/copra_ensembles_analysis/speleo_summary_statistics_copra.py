# -*- coding: utf-8 -*-
#
# Copyright (C) 2008-2014 Jonathan F. Donges
# Author: Jonathan F. Donges <donges@pik-potsdam.de>
# URL: <http://www.pik-potsdam.de/members/donges/software>

"""
Computes and plots a summary statistics counting the relative number of events
in paleoclimate variability occuring in several records at the same time using
COPRA ensembles.

This script provides analyses for this publication:

J.F. Donges, R.V. Donner, N. Marwan, S.F.M. Breitenbach, K. Rehfeld,
and J. Kurths,
Nonlinear regime shifts in Holocene Asian monsoon variability:
Potential impacts on cultural change and migratory patterns,
Climate of the Past 11, 709-741 (2015),
DOI: 10.5194/cp-11-709-2015

For the COPRA algorithm refer to:

S.F.M. Breitenbach  et al. (2012),
COnstructing Proxy Records from Age models (COPRA),
Climate of the Past 8, 1765–1779, doi: 10.5194/cp-8-1765-2012.
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

#  Import progress bar for easy progress bar handling
import progressbar


#
#  Settings
#

#  Settings for statistics

#  Quantiles for statistical testing

# Used for surrogates of measure and summary statistics
Q_QUANTILE = 0.10
# Used for defining confidence bounds from COPRA ensemble
Q_QUANTILE_COPRA = 0.10

#  Settings for plotting

#  List of letters from "A" to "Z"
figure_labels = map(chr, range(65, 91))

#  Toggle plotting of confidence band
SHOW_CONFIDENCE_BAND = True

#  Set labels
AGE_LABEL = "Age (ka B.P.)"
DATA_LABEL = r"$\Delta \delta^{18}$O" #  (‰ VPDB)

#  Settings for plotting additional information

# See Mayewski et al., Quat. Sci. Rev. (2004)
BOND_EVENTS = np.array([10.3, 9.4, 8.2, 5.9, 4.35, 2.8, 1.2, 0.5])
RCC_EPISODES = np.array([[8., 9.], [6., 5.], [3.8, 4.2], [2.5, 3.5], [1., 1.2],
                         [0.15, 0.6]])

#  Settings for computation of summary statistics

INTERP_TIME_RES = 10. # years [a]
SIGN_WINDOW_SIZE = 100. # years [a]
# Number of surrogates for Monte Carlo significance testing of summary
# statistics
N_SURROGATES_SUMM_STAT = 1000

#
#  List of data files
#

#  List files in the desired plotting order
DATA_FILE_NAMES = ["results_copra_ensemble_DIM_3_TAU_25_M_100_DETREND_True_DIMARSHIM.pickle",
                   "results_copra_ensemble_DIM_3_TAU_54_M_100_DETREND_True_QUNF.pickle",
                   "results_copra_ensemble_DIM_3_TAU_14_M_100_DETREND_True_HOTI.pickle",
                   "results_copra_ensemble_DIM_3_TAU_36_M_100_DETREND_True_MAWMLUH.pickle",
                   "results_copra_ensemble_DIM_3_TAU_22_M_100_DETREND_True_TIANMEN.pickle",
                   "results_copra_ensemble_DIM_3_TAU_46_M_100_DETREND_True_DONGGE.pickle",
                   "results_copra_ensemble_DIM_3_TAU_48_M_100_DETREND_True_LIANHUA.pickle",
                   "results_copra_ensemble_DIM_3_TAU_15_M_100_DETREND_True_HESHANG.pickle",
                   "results_copra_ensemble_DIM_3_TAU_18_M_100_DETREND_True_JIUXIAN.pickle",
                   "results_copra_ensemble_DIM_3_TAU_33_M_100_DETREND_True_LIANG-LUAR1.pickle",
                   "results_copra_ensemble_DIM_3_TAU_33_M_100_DETREND_True_LIANG-LUAR2.pickle",
                   "results_copra_ensemble_DIM_3_TAU_33_M_100_DETREND_True_LIANG-LUAR3.pickle"]

n_time_series = len(DATA_FILE_NAMES)


#
#  Load results for plotting
#

#  Initialize dicts
storage = {}
FILENAME = {}
NAME = {}
N_ENSEMBLE = {}
symbols = {}
time = {}
values = {}
step_sequence = {}
results = {}
surrogate_results = {}

#  Load from files
for i in xrange(n_time_series):
    print "Loading results from", DATA_FILE_NAMES[i]

    file = open(DATA_FILE_NAMES[i], 'r')
    storage[i] = cPickle.load(file)
    file.close()

    #  Store parameters
    FILENAME[i] = storage[i]["FILENAME"]
    NAME[i] = storage[i]["NAME"]
    N_ENSEMBLE[i] = storage[i]["N_ENSEMBLE"]

    #  Store symbols
    #  !!! Assume that the same symbols are used in all files!
    symbols = storage[i]["symbols"]

    #  Store raw input data
    time[i] = storage[i]["time"]
    values[i] = storage[i]["values"]

    #  Store axes
    step_sequence[i] = storage[i]["step_sequence"]

    #  Store results
    results[i] = storage[i]["results"]
    surrogate_results[i] = storage[i]["surrogate_results"]

    #  Derived properties
    print "Number of realizations for significance test:", N_ENSEMBLE[i]

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

#  Init dicts
q_elements = {}
median = {}
q_left = {}
q_right = {}
surr_q_left = {}
surr_q_right = {}

for i in xrange(n_time_series):
    #  Number of elements right and left of the quantile
    #  Divide by two, such that the total probability of false rejection of
    #  null hypothesis
    #  is equal to Q_QUANTILE_COPRA (two sided test).
    q_elements[i] = int(Q_QUANTILE_COPRA / 2 * N_ENSEMBLE[i])

    #  Get statistics for results
    median[i] = {}
    q_left[i] = {}
    q_right[i] = {}

    for measure in symbols.keys():
        results[i][measure].sort(axis=0)

        #  Get median
        median[i][measure] = np.median(results[i][measure], axis=0)

        #  Get q-quantiles
        q_left[i][measure] = results[i][measure][q_elements[i] - 1, :]
        q_right[i][measure] = results[i][measure][-q_elements[i], :]


    #
    #  Calculate q-quantiles for surrogates
    #

    q_elements[i] = int(Q_QUANTILE / 2 * N_ENSEMBLE[i])

    surr_q_left[i] = {}
    surr_q_right[i] = {}

    for measure in symbols.keys():
        #  Get quantiles from surrogates from all time series realizations for now
        surrogate_results[i][measure] = surrogate_results[i][measure].flatten()

        #  Sort
        surrogate_results[i][measure].sort()

        #  Get q-quantiles
        surr_q_left[i][measure] = surrogate_results[i][measure][q_elements[i] - 1]
        surr_q_right[i][measure] = surrogate_results[i][measure][-q_elements[i]]


#
#  Prepare time axis
#

min_age_list = []
max_age_list = []

for i in xrange(n_time_series):
    #  Get minimum and maximum age in time series
    min_age_list.append(time[i].min())
    max_age_list.append(time[i].max())

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
    median_interp = np.zeros((n_time_series, len(time_interp)))
    q_left_interp = np.zeros((n_time_series, len(time_interp)))
    q_right_interp = np.zeros((n_time_series, len(time_interp)))
    is_significant = np.zeros((n_time_series, len(time_interp)), dtype=bool)
    event_lengths = np.zeros((n_time_series, len(time_interp)), dtype="int32")

    for i in xrange(n_time_series):
        #  Interpolate result median and quantiles piecewise linearly
        median_interp[i,:] = np.interp(time_interp, time[i][step_sequence[i]],
                                       median[i][measure], left=np.nan,
                                       right=np.nan)
        q_left_interp[i,:] = np.interp(time_interp, time[i][step_sequence[i]],
                                       q_left[i][measure], left=np.nan,
                                       right=np.nan)
        q_right_interp[i,:] = np.interp(time_interp, time[i][step_sequence[i]],
                                        q_right[i][measure], left=np.nan,
                                        right=np.nan)

        #  Determine at which time steps the result quantiles are outside the surrogate confidence band indicating events/transitions
        is_significant[i,:] = np.logical_or(q_left_interp[i,:] < surr_q_left[i][measure], q_right_interp[i,:] > surr_q_right[i][measure])

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
    n_records_at_t = n_time_series - np.isnan(median_interp).sum(axis=0)

    #  Compute summary statistics
    (max_records_around_t, significant_relative[measure]) = summary_statistics(is_significant, time_interp, n_records_at_t, n_window)

    #  Store results
    np.savetxt("time_axis.txt", time_interp)
    np.savetxt("max_records_around_t.txt", max_records_around_t)
    np.savetxt("summary_stat_" + measure + ".txt",
               significant_relative[measure])

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
    event_surrogates = np.zeros((n_time_series, N_SURROGATES_SUMM_STAT,
                                len(time_interp)), dtype=bool)

    for i in xrange(n_time_series):
        #  Get indices corresponding to time boundaries of present record
        min_time = time[i][step_sequence[i]][0]
        max_time = time[i][step_sequence[i]][-1]

        min_time_index = min(range(len(time_interp)),
                             key=lambda i: abs(time_interp[i] - min_time))
        max_time_index = min(range(len(time_interp)),
                             key=lambda i: abs(time_interp[i] - max_time))

        #  Repeat until desired number of valid surrogates is generated
        surr_count = 0
        while surr_count < N_SURROGATES_SUMM_STAT:
            #  Create surrogate

            #  Loop over event lengths
            for j in xrange(1, len(time_interp)):
                #  Loop over all events of a given length j
                for k in xrange(event_length_dist[measure][i,j]):
                    #  Randomly draw index where event of length j begins and
                    #  make sure that event cannot exceed time span of
                    #  original record
                    index = int(np.random.rand() * (max_time_index - min_time_index - j) + min_time_index)

                    event_surrogates[i, surr_count, index:index + j] = True

            #  Count surrogate
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
        summary_stat[j,:] = summary_statistics(event_surrogates[:,j,:],
                                               time_interp, n_records_at_t,
                                               n_window)[1]

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

for time in time:
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
pylab.fill_between(time_interp / 1000., 0, max_records_around_t,
                   color="#F1B08E")

pylab.ylabel(r"$N_r$", fontsize=16)

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

    #  Prepare lower bound for fill between to highlight episodes with
    #  significant excursions from "normal" climate state
    is_significant = q_sum_stat_upper[measure] < significant_relative[measure]

    #  Plot results
    pylab.plot(time_interp / 1000., significant_relative[measure],
               colors_line[j - 2])
    #  Fill insignificant episodes
    pylab.fill_between(time_interp / 1000., 0, significant_relative[measure],
                       color=colors_fill[j - 2])
    #  Fill significant episodes
    pylab.fill_between(time_interp / 1000., q_sum_stat_upper[measure],
                       significant_relative[measure], where=is_significant,
                       color=colors_line[j - 2])

    #  Plot confidence band for summary statistics
    #pylab.fill_between(time_interp / 1000., q_sum_stat_lower[measure], q_sum_stat_upper[measure], color="0.75")
    #pylab.plot(time_interp / 1000., q_sum_stat_lower[measure], "r-")
    pylab.plot(time_interp / 1000., q_sum_stat_upper[measure], "-",
               color="grey")

    pylab.ylabel("$n_" + symbols[measure][1:-1] + "$", fontsize=16)

    #  Plot additional information
    #  Plot RCC episodes
    for k in xrange(RCC_EPISODES.shape[0]):
        pylab.axvspan(RCC_EPISODES[k,0], RCC_EPISODES[k,1], fill=True,
                      edgecolor="0.8", color="0.8", zorder=0)
    #  Plot Bond events
    pylab.scatter(BOND_EVENTS, 0.9 * np.ones(len(BOND_EVENTS)), c='k',
                  marker='o', s=80)

    pylab.xlim(0, 10)
    pylab.ylim(0, 1)

    #  Add figure label
    ax.annotate(figure_labels[j - 1], xy=(-0.14, 1.1),
                xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top',
                fontname='Helvetica', fontsize=14, fontweight='bold')

    #  Render time labels invisible
    #ax.axes.get_xaxis().set_visible(False)

    #  Enable minor ticks
    ax.minorticks_on()

    #  Render time tick labels invisible
    ax.set_xticklabels([])

    #  Set number of ticks and tick labels
    #ax.locator_params(axis="x", tight=True, nbins=9)
    ax.locator_params(axis="y", tight=True, nbins=5)

#  Enable minor ticks
ax.minorticks_on()

pylab.xlabel(AGE_LABEL)

#  Save figure
pylab.savefig("summary_statistics_copra.pdf")

#  Show all figures
pylab.show()
