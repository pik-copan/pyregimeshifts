# -*- coding: utf-8 -*-
#
# Copyright (C) 2008-2014 Jonathan F. Donges
# Author: Jonathan F. Donges <donges@pik-potsdam.de>
# URL: <http://www.pik-potsdam.de/members/donges/software>

"""
Computes auto-correlation function for irregularly sampled time series.

Uses the method proposed in:

Rehfeld, K., Marwan, N., Heitzig, J., & Kurths, J. (2011). Comparison of correlation analysis techniques for irregularly sampled time series. Nonlinear Processes in Geophysics, 18(3), 389-404.

This script provides analyses for this publication:

J.F. Donges, R.V. Donner, N. Marwan, S.F.M. Breitenbach, K. Rehfeld, and J. Kurths,
Nonlinear regime shifts in Holocene Asian monsoon variability: Potential impacts on cultural change and migratory patterns,
Climate of the Past 11, 709-741 (2015),
DOI: 10.5194/cp-11-709-2015
"""

#
#  Imports
#

import sys

import numpy as np

import pylab
import progressbar

#
#  Settings
#

#  Filename
FILENAME_X = "../../data/raw_proxy_data/Dongge_DA.dat"

#  Resolution of cross-correlation (units of time)
DELTA_LAG = 10 #  Measured in years here

#  Maximum lag index
MAX_LAG_INDEX = 100

#  Toggle detrending
DETRENDING = True
DETRENDING_WINDOW_SIZE = 1000. # Measured in years here


#
#  Functions
#

def detrend_time_series(time, data, window_size):
    #  Get length of data array
    n = data.shape[0]
    #  Initialize a local copy of data array
    detrended_data = np.empty(n)

    #  Detrend data
    for j in xrange(n):
        #  Get lower and upper bound of window in time domain
        lower_bound = time[j] - window_size / 2.
        upper_bound = time[j] + window_size / 2.

        #  Get time indices lying within the window
        window_indices = np.logical_and(time >= lower_bound, time <= upper_bound)

        #  Substract window mean from data point in the center
        detrended_data[j] = data[j] - data[window_indices].mean()

    return detrended_data

def gaussian(x, std):
    """
    Returns value of gaussian distribution at x with 0 mean
    and standard deviation std.
    """
    return 1 / np.sqrt(2 * np.pi * std) * np.exp(-np.abs(x ** 2) / (2 * std**2) )

def kernel_auto_correlation_est(x, time_diff, kernel_func, kernel_param,
                                 delta_lag, max_lag_index):
    """
    Estimates auto correlation using a kernel function.
    """
    #  Normalize time series
    x -= x.mean()
    x /= x.std()

    #  Initialize discrete auto-correlation function
    auto_correlation = np.zeros(max_lag_index + 1)

    #  Loop over all positive lags and zero lag
    for k in xrange(max_lag_index + 1):
        #  Calculate b matrix
        b = kernel_func(k * delta_lag - time_diff, kernel_param)

        #  Calculate nominator
        nominator = np.dot(x, np.dot(b, x.transpose()))
        #  Calculate denominator
        denominator = b.sum()

        #  Calculate auto-correlation
        auto_correlation[k] = nominator / denominator

    lag_times = delta_lag * np.arange(max_lag_index + 1)

    return (lag_times, auto_correlation)

#
#  Main script
#

#  Load record x
data_x = np.loadtxt(FILENAME_X, unpack=False, usecols=(0,1,), comments="#")
#data_x = np.fromfile(FILENAME_X, sep=" ")
time_x = data_x[:,0]
x = data_x[:,1]

#  Detrending of time series using moving window averages
if DETRENDING:
    x = detrend_time_series(time_x, x, DETRENDING_WINDOW_SIZE)

#  Get length of records
N_x = len(time_x)

#  Get recommended standard deviation of gaussian Kernel (Kira Rehfeld's
#  NPG paper)
sigma = 0.25 * np.diff(time_x).mean()

print "Length of record x:", N_x
print "Mean sampling time x:", np.diff(time_x).mean()
print "Recommended standard deviation of gaussian Kernel:", sigma

#  Calculate matrix of time differences
time_diff = np.zeros((N_x, N_x))

for i in xrange(N_x):
    for j in xrange(N_x):
        time_diff[i,j] = time_x[i] - time_x[j]

#  Estimate auto-correlation function
(lag_times, auto_correlation) = kernel_auto_correlation_est(x=x.copy(), time_diff=time_diff, kernel_func=gaussian, kernel_param=sigma, delta_lag=DELTA_LAG, max_lag_index=MAX_LAG_INDEX)


#
#  Save results
#

results = np.zeros((MAX_LAG_INDEX + 1, 2))

results[:,0] = lag_times
results[:,1] = auto_correlation

np.savetxt("kernel_acf_dongge.txt", results)

#
#  Plot results
#

#  Set plotting parameters (for Clim. Past paper)
params = { 'figure.figsize': (6.,6.),
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

#  Plot time series
pylab.figure(1)

pylab.plot(time_x, x)
pylab.xlabel("Age (y B.P.)")
pylab.ylabel("Normalized values")

pylab.figure(2)

pylab.plot(lag_times, auto_correlation, "k")

pylab.axhline(y=1 / np.e, color="red")
pylab.xlabel("Time delay [y]")
pylab.ylabel("ACF")

pylab.ylim(-0.5,1)

pylab.savefig("auto_corr_irregular.pdf")

pylab.show()
