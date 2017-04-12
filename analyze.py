#!/usr/bin/python
# A script to analyze synchrony in neural networks
from __future__ import print_function


from numpy import *
from pylab import *



# ENTRY CONDITIONS:
# - `means` is a list of the mean firing rate over time
# EXIT CONDITIONS:
# - Draws the original timeseries, an FFT, and the autocorrelation
def draw_figure(means):
    subplot(3, 1, 1)
    title("Means")
    plot(means)

    subplot(3, 1, 2)
    title("FFT")
    plot(abs(fft(means)))

    subplot(3, 1, 3)
    title("Autocorrelation")
    plot(_autocorrelation(_smooth(means)))
    show()

# ENTRY CONDITIONS:
# - `means` is a list of the mean firing rate over time
# EXIT CONDITIONS:
# - Prints out a number of analysis measures, most of which are only
#   experimental at the moment
def print_analysis_measures(means):
    # A higher number here should signify high synchrony
    print("Ratio of the amplitude of the second peak to the first:", analyze_amplitude_ratio(means))

    # A higher number here should signify high synchrony
    print ("Ratio of the amplitude of the second peak to the first in autocorrelation:", analyze_autocor(means))

    #High individual numbers
    fft_spectrum = abs(fft(means))
    fft_peaks = _peakdet(fft_spectrum)
    print ("Ratio of the amplitude of each fft peak to the mean:")
    for i in range(0, len(fft_peaks)):
        print ("\tPeak at %i: %f" % (fft_peaks[i], fft_spectrum[fft_peaks[i]]/mean(fft_spectrum)))
    
# ENTRY CONDITIONS:
# - `ts` is a timeseries.  Usually this is the list of the mean firing
#   rates over time.
# EXIT CONDITIONS:
# - Returns the amplitude of the second peak normalized by the maximum
#   amplitude of the curve.  If only one peak is present (or no
#   peaks), 0 is returned.

def analyze_amplitude_ratio(ts):
    peaks = _peakdet(ts)
    max_amplitude = max(ts) - min(ts)
    if len(peaks) < 3:
        return -1
    amp = ts[peaks[2]]-ts[peaks[1]]
    return amp/mean(ts)

# This is just a shortcut function for:
#     analyze_amplitude_ratio(_autocorrelation(means))

# The main purpose of this is to keep _autocorrelation a private
# method and to reinforce that this is one of our core methods of
# analysis.  It will only analyze the first half of the
# autocorrelation spectrum, since this is the only part that is really
# relevant.
def analyze_autocor(means):
    return analyze_amplitude_ratio(_autocorrelation(_smooth(means))[0:int(len(means)/2)])

# ENTRY CONDITIONS:
# - *means* is a list of the mean firing rate over time
# EXIT CONDITIONS
# - Returns the autocorrelation of `means`.  The autocorrelation is
#   the covariance of the data points with the points `delta` in front
#   of them as an array indexed by `delta`.
# def _autocorrelation(means):
#     l = len(means)
#     cov = []
#     for i in range(0,l):
#         array1 = means[:l-i]
#         array2 = means[i:]
#         mean1 = average(array1)
#         mean2 = average(array2)
#         covar = sum(multiply((array1-mean1), (array2-mean2)))/(l-i)
#         cov.append(covar/l)
#     return cov
def _autocorrelation(means):
    from numpy.core import multiarray
    result = multiarray.correlate(means, means, mode=2)
    return result[result.size/2:]

# ENTRY CONDITIONS:
# - `ts` is a timeseries to smooth
# - `window` is the window by which to smooth it
# EXIT CONDITIONS:
# - Returns a timeseries that is a smoother version of `ts`
def _smooth(ts, window=10):
    finallist = []
    for i in range(0, len(ts)-window):
        finallist.append(average(ts[i:i+window]))
    finallist.extend(ts[len(ts)-window:len(ts)-1])
    return finallist

# ENTRY CONDITIONS:
# - *ts* is a list of float numbers
# - Optionally, *threshold_ratio* is a floating point threshold for
#   the minimum amplitude of a peak, compared to the highest value in
#   the spectrum minus the lowest value
# EXIT CONDITIONS:
# - Returns a list of indices of *ts* that represent local maxima and
#   minima.  Odd indices are minima, even indices are maxima.

def _peakdet(ts, threshold_ratio=.1):
    """A peak detection algorithm.  This works using the knowledge
    that a maximum must come in between two minima, and vice versa.
    When it finds a suspected maximum value, it saves the index and
    continues scanning the timeseries until it finds either a greater
    value, or a value that is *threshold* less than it.  Upon finding
    the latter, the previous greatest value is saved as a peak, and
    the algorithm proceeds to look for a minimum in the same manner."""
    THRESH = threshold_ratio * (max(ts)-min(ts))
    maxima = []
    minima = []
    extrema = []
    looking_for_maximum = True
    last = 0
    for i in range(1, len(ts)):
        if looking_for_maximum:
            if ts[i] > ts[last]:
                last = i
            elif ts[i] + THRESH < ts[last]:
                maxima.append(last)
                extrema.append(last)
                looking_for_maximum = False
        else: #looking for minimum
            if ts[i] < ts[last]:
                last = i
            elif ts[i] - THRESH > ts[last]:
                minima.append(last)
                extrema.append(last)
                looking_for_maximum = True
        
    return extrema
        
if __name__ == "__main__":
    series = multiply(sin(pow(linspace(1, 100, 1000), .9)), exp(linspace(1, -10, 1000)))
    draw_figure(series)
    print_analysis_measures(series)
    show()
    #draw_figure(exp(linspace(1, 100, 1000)))
    #show()
