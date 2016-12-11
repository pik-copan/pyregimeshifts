# -*- coding: utf-8 -*-
#
# Copyright (C) 2008-2014 Jonathan F. Donges
# Author: Jonathan F. Donges <donges@pik-potsdam.de>
# URL: <http://www.pik-potsdam.de/members/donges/software>

"""
Performs recurrence analysis of paleoclimate proxy records.

This script provides analyses for this publication:

J.F. Donges, R.V. Donner, N. Marwan, S.F.M. Breitenbach, K. Rehfeld, and J. Kurths,
Nonlinear regime shifts in Holocene Asian monsoon variability: Potential impacts on cultural change and migratory patterns,
Climate of the Past 11, 709-741 (2015),
DOI: 10.5194/cp-11-709-2015
"""

#
#  Imports
#

#  Import cPickle for loading and saving data
import cPickle

#  Import np for fast numerics
import numpy as np

#  Import progress bar for easy progress bar handling
import progressbar

#  Import class for recurrence network analysis
from pyunicorn.timeseries import RecurrenceNetwork

#
#  Settings
#

#  Name of data directory
DATA_DIR = "../../data/raw_proxy_data/"

#  List of data FILENAMES
FILENAMES = ["Dimarshim_D1.dat", "Qunf_Q5_orig.dat", "Hoti.dat",
             "Mawmluh.dat", "Tianmen_TM18_older.dat", "Dongge_DA.dat",
             "Lianhua_d18O_d13C.dat", "Heshang_HS4.dat", "Jiuxian.dat",
             "Liang-Luar.dat"]

#  Names of proxy records / caves
NAMES = ["Dimarshim", "Qunf", "Hoti", "Mawmluh", "Tianmen", "Dongge",
         "Lianhua", "Heshang", "Jiuxian", "Liang-Luar"]


#  Specify symbol used for commenting in data file
COMMENT_SYMBOL = "%"

#  Settings for the time dependent recurrence plot

#  Window length [a] / [ka]
T_TIME = [750., 750., 750., 750., 750., 750., 750., 750., 750. ,750.]

#  Step size [a] / [ka]
DELTA_TIME = 50.

#  Settings for the embedding
DIM = 3
TAU = 2 # Only used if ADAPT_DELAY == False

ADAPT_DELAY = True #  If true, the delay in units of data points is estimated to match the given DELAY_TIME

#  Explicitly set delay times for time-delay embedding
DELAY_TIMES = [100., 216., 57., 146., 90., 185., 193., 60., 73., 135.] # in years [a]

#  Settings for the recurrence plot
METRIC = "supremum" # metric for recurrence definition
RR = 0.05 # prescribed recurrence rate

#  Settings for significance testing

#  Ensemble size
N_ENSEMBLE = 1000
#  Choose whether whole embedded state vectors or the scalar time series should be shuffled (Different null-hypothesis!)
SHUFFLE_EMBEDDED = True

#  Settings for detrending
DETREND = True
DETRENDING_WINDOW_SIZE = 1000. # measured in [a] / [ka]

#
#  Functions
#

def detrend_time_series(data, window_size):
    """
    """
    #  Get length of data array
    n = data.shape[0]
    #  Initialize a local copy of data array
    detrended_data = np.empty(n)

    #  Detrend data
    for j in xrange(n):
        #  Get distance of sample from boundaries of time series
        dist = min(j, n - 1 - j)

        if window_size / 2 > dist:
            half_size = dist
        else:
            half_size = window_size / 2

        detrended_data[j] = data[j] - data[j - half_size:j + half_size + 1].mean()

    return detrended_data

def autocorrelation(data, lag):
    """Return autocorrelation of data at specified lag."""
    return np.corrcoef(data[lag:], data[:-lag])[0,1]

#
#  The main script
#

print "Recurrence network analysis of paleoclimate records"
print "---------------------------------------------------"

#
#  Import data
#

time_list = []
data_list = []
sampling_time_list = []
delay_list = []

#  Get number of time series
n_time_series = len(FILENAMES)

#  Load data
for i in xrange(n_time_series):
    time, data = np.loadtxt(DATA_DIR + FILENAMES[i], comments=COMMENT_SYMBOL,
                            unpack=True, usecols=(0,1,))
    average_sampling_time = np.diff(time).mean()

    time_list.append(time)
    if DETREND:
        #  Detrend data!
        detrended_data = detrend_time_series(data=data,
                window_size=DETRENDING_WINDOW_SIZE / average_sampling_time)
        data_list.append(detrended_data)

    else:
        data_list.append(data)

    #  Get average sampling times
    sampling_time_list.append(average_sampling_time)

    #  Get delay time
    delay_list.append(int(DELAY_TIMES[i] / average_sampling_time))

    #  Temporaray: Get length of time series
    n = len(time)


#
#  Print some statistics
#

print "Average sampling time:"

for i in xrange(n_time_series):
    print FILENAMES[i], ": (", np.diff(time_list[i]).mean(), "pm", np.diff(time_list[i]).std(), ") ka"

#
#  Analyze time dependent recurrence networks by moving a window over
#  the time series
#

#  Initialize list of window mid-points used for estimating time scale of windowed measures
step_sequence = []

#  Create dictionary of symbols for each windowed measure to be calculated
symbols = {"Average path length": "$\mathcal{L}$",
           "Transitivity": "$\mathcal{T}$"}

#symbols = {"Average path length": "$\mathcal{L}$",
#           "n.s.i. average path length": "$\mathcal{L}^*$",
#           "Clustering": "$\mathcal{C}$",
#           "n.s.i. clustering": "$\mathcal{C}^*$"}

#symbols = {"Determinism": "$DET$",
#           "Laminarity": "$LAM$",
#           "Mean diagonal line length": "$L_{mean}$",
#           "Trapping time": "$TT$",
#           "Diagonal line entropy": "$ENTR$",
#           "Autocorrelation": "$ACF(1)$",
#           "Mean": "Mean",
#           "Standard deviation": "STD"}

#  Initialize dictionaries
results = {}
surrogate_results = {}

for measure in symbols.keys():
    results[measure] = []
    surrogate_results[measure] = []

#  Run analysis for each time series separately
for i in xrange(n_time_series):
    print "Analyzing original data from", FILENAMES[i]

    #  Get time and data arrays
    time = time_list[i]
    data = data_list[i]
    sampling_time = sampling_time_list[i]

    #  Set delay
    if ADAPT_DELAY:
        TAU = delay_list[i]

    #  Get window and step size in units of samples
    T = int(T_TIME[i] / sampling_time)
    delta = int(DELTA_TIME / sampling_time)

    #  Get length of time series
    t_max = len(time)

    #  Get required time series length before embedding to achive window length T in the recurrence plot
    T_embedded = T + (DIM - 1) * TAU

    #  Get number of steps
    t_steps = int((t_max - T_embedded) / float(delta) + 1)

    print "Length of record:", t_max
    print "Size of moving window:", T
    print "Step size:", delta
    print "Number of steps for moving window:", t_steps
    print "Embedding dimension:", DIM
    print "Embedding delay:", TAU
    print "Prescribed link density / recurrence rate:", RR

    #  Initializations
    local_step_sequence = np.empty((t_steps), dtype=int)
    local_result = {}

    for measure in symbols.keys():
        local_result[measure] = np.empty(t_steps)

    #  Initialize progress bar
    progress = progressbar.ProgressBar().start()

    #  Loop over moving windows
    for j in xrange(t_steps):
        #  Get time series section for current window
        time_series = data[j * delta:j * delta + T_embedded]
        local_step_sequence[j] = j * delta + T_embedded / 2

        #  Prepare recurrence network from original data
        rec_net = RecurrenceNetwork(time_series.flatten(), dim=DIM, tau=TAU,
                                    metric=METRIC, normalize=False,
                                    silence_level=2, recurrence_rate=RR)

        #  Calculations for original recurrence network
        local_result["Average path length"][j] = rec_net.average_path_length()
        local_result["Transitivity"][j] = rec_net.transitivity()

        #local_result["Assortativity"][j] = rec_net.assortativity()
        #local_result["Diameter"][j] = rec_net.diameter()

        #  Calculate RQA measures
        #local_result["Determinism"][j] = rec_net.determinism()
        #local_result["Laminarity"][j] = rec_net.laminarity()
        #local_result["Mean diagonal line length"][j] = rec_net.average_diaglength()
        #local_result["Trapping time"][j] = rec_net.trapping_time()
        #local_result["Diagonal line entropy"][j] = rec_net.diag_entropy()
        #local_result["Autocorrelation"][j] = autocorrelation(time_series, lag=1)
        #local_result["Mean"][j] = time_series.mean()
        #local_result["Standard deviation"][j] = time_series.std()

        #  Update progress bar every step
        progress.update(int(100 * j / float(t_steps)))

    #  Terminate progress bar
    progress.finish()

    #  Store window mid-point
    step_sequence.append(local_step_sequence)

    #  Store results
    for measure in symbols.keys():
        results[measure].append(local_result[measure])

    #
    #  Calculate significance levels for network measures
    #

    print "Calculating significance levels based on", N_ENSEMBLE, "surrogates..."

    #  Initialize progress bar
    progress = progressbar.ProgressBar().start()

    #  Create a copy of data for generating surrogates from
    surrogate_data = data.copy()

    if SHUFFLE_EMBEDDED:
        #  Get embedding of full time series
        surrogate_embedding = rec_net.embed_time_series(surrogate_data,
                                                        DIM, TAU)

    #  Prepare stuff
    local_surrogate_result = {}

    for measure in symbols.keys():
        local_surrogate_result[measure] = np.empty(N_ENSEMBLE)

    for j in xrange(N_ENSEMBLE):
        if SHUFFLE_EMBEDDED:
            #  Shuffle embedded time series along time axis, that is, whole
            #  embedded state vectors are shuffled around.
            permuted_indices = np.random.permutation(surrogate_embedding.shape[0])

            #  Use the first T state vectors from the shuffled and embedded
            #  time series as a surrogate for one window
            surrogate_series = surrogate_embedding[permuted_indices[:T],:]

            #  Prepare recurrence network from surrogate data for shuffled
            #  embedded time series
            rec_net = RecurrenceNetwork(surrogate_series.copy(),
                                        metric=METRIC, normalize=False,
                                        silence_level=2, recurrence_rate=RR)
        else:
            #  Shuffle dust time series
            permuted_indices = np.random.permutation(surrogate_data.shape[0])

            #  Use the first T_embedded states from the shuffled dust time series as a surrogate for one window
            surrogate_series = surrogate_data[permuted_indices[:T_embedded]]

            #  Prepare recurrence network from surrogate data for shuffled time series
            rec_net = RecurrenceNetwork(surrogate_series.copy(), dim=DIM,
                                        tau=TAU, metric=METRIC,
                                        normalize=False, silence_level=2,
                                        recurrence_rate=RR)

        #  Calculate measures for surrogate network
        local_surrogate_result["Average path length"][j] = rec_net.average_path_length()
        local_surrogate_result["Transitivity"][j] = rec_net.transitivity()

        #local_surrogate_result["Assortativity"][j] = rec_net.assortativity()
        #local_surrogate_result["Diameter"][j] = rec_net.diameter()

        #  Calculate RQA measures
        #local_surrogate_result["Determinism"][j] = rec_net.determinism()
        #local_surrogate_result["Laminarity"][j] = rec_net.laminarity()
        #local_surrogate_result["Mean diagonal line length"][j] = rec_net.average_diaglength()
        #local_surrogate_result["Trapping time"][j] = rec_net.trapping_time()
        #local_surrogate_result["Diagonal line entropy"][j] = rec_net.diag_entropy()
        #local_surrogate_result["Autocorrelation"][j] = autocorrelation(data, lag=1)
        #local_surrogate_result["Mean"][j] = data.mean()
        #local_surrogate_result["Standard deviation"][j] = data.std()

        #  Update progress bar every step
        progress.update(int(100 * j / float(N_ENSEMBLE)))

    #  Store results
    for measure in symbols.keys():
        surrogate_results[measure].append(local_surrogate_result[measure])

    #  Terminate progress bar
    progress.finish()

#
#  Save results
#

print "Saving results..."

#  Initialize storage dictionary
storage = {}

#  Store parameters
storage["FILENAMES"] = FILENAMES
storage["NAMES"] = NAMES

storage["T_TIME"] = T_TIME
storage["DELTA_TIME"] = DELTA_TIME
storage["DETRENDING_WINDOW_SIZE"] = DETRENDING_WINDOW_SIZE

storage["DIM"] = DIM
storage["TAU"] = TAU
storage["ADAPT_DELAY"] = ADAPT_DELAY
storage["DELAY_TIMES"] = DELAY_TIMES

storage["METRIC"] = METRIC
storage["RR"] = RR

storage["N_ENSEMBLE"] = N_ENSEMBLE
storage["SHUFFLE_EMBEDDED"] = SHUFFLE_EMBEDDED

#  Store symbols
storage["symbols"] = symbols

#  Store raw input data
storage["time_list"] = time_list
storage["data_list"] = data_list

#  Store axes
storage["step_sequence"] = step_sequence

#  Store results
storage["results"] = results
storage["surrogate_results"] = surrogate_results

#  Save to file
filename = "results_speleo_comparison_W_" + str(T_TIME[0]) + "y_M_" + str(N_ENSEMBLE) + "_DETREND_" + str(DETREND) + ".pickle"
file = open(filename, 'w')
cPickle.dump(storage, file)
file.close()
