"""
Compute intrinsic neural timescale (intsc) from time-series data.

Gao, R., van den Brink, R. L., Pfeffer, T., & Voytek, B. (2020). Neuronal
timescales are functionally dynamic and shaped by cortical microarchitecture.
eLife 2020;9:e61277

Example usage

imgfile=~/data/NDARAM848GTE_restingstate.txt
outfile=~/data/NDARAM848GTE_restingstate_intsc
srate=500

python intsc_gao.py -i $imgfile -o $outfile --srate 500
"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd
from fooof import FOOOF, FOOOFGroup
from optparse import OptionParser


# function to parse input arguments
def parse_args():
    """
    Parse arguments.
    """

    parser=OptionParser()

    parser.add_option('-i',"--input", \
                      dest='imgfile', \
                      help="Filename of the input EEG dataset", \
                      default=None)
    parser.add_option('-o',"--output", \
                      dest='outname', \
                      help="Output filename prefix", \
                      default=None)
    parser.add_option('',"--srate", \
                      dest='srate',\
                      help="sampling rate in Hz", \
                      default=None)

    (options,args) = parser.parse_args()

    return(options)


def convert_knee_val(knee, exponent=2.):
    """
    Convert knee parameter to frequency and time-constant value.
    Can operate on array or float.
    Default exponent value of 2 means take the square-root, but simulation shows
    taking the exp-th root returns a more accurate drop-off frequency estimate
    when the PSD is actually Lorentzian.
    """
    knee_freq = knee**(1./exponent)
    knee_tau = 1./(2*np.pi*knee_freq)
    return knee_freq, knee_tau


# function to load data
def load_data(datafile):
    """
    Load data from a file
    """

    # read in data
    data = pd.read_csv(datafile, delimiter = '\t', index_col = 0)

    # check if columns have NaNs
    colnames = data.columns
    nrows = data.shape[0]
    ncols = data.shape[1]
    mask = np.full(shape = ncols, fill_value = True, dtype=bool)
    for col_idx, column in enumerate(colnames):
        n_nan = np.sum(np.isnan(data.loc[:,column]))
        if n_nan == nrows:
            mask[col_idx] = False

    # grab only columns with no NaNs
    data = data.iloc[:,mask]
    return(data)


# main function for computing intrinsic neural timescale
def compute_intsc(ts, numtps, srate):
    """
    Compute intrinsic neural timescale
    """

    # figure out tr in seconds
    tr = 1/srate

    # grab electrode names
    colnames2use = ts.columns
    colnames2use = colnames2use.to_list()

    # normalize (divide all elements of the timeseries by the first value in the timeseries)
    ts = np.array(ts)
    for col_index, ts2use in enumerate(ts.T):
        ts[:,col_index] = np.divide(ts2use,ts2use[0])

    # compute spectrum via FFT
    f_axis = np.fft.fftfreq(numtps, tr)[:int(np.floor(numtps/2))]
    psds = (np.abs(sp.fft(ts, axis=0))**2)[:len(f_axis)]

    # fit FOOOF & get knee parameter & convert to timescale
    fooof = FOOOFGroup(aperiodic_mode = 'knee', max_n_peaks = 0, verbose = False)
    fooof.fit(freqs = f_axis, power_spectra = psds.T, freq_range = (2,200))
    fit_knee = fooof.get_params('aperiodic_params', 'knee')
    fit_exp = fooof.get_params('aperiodic_params', 'exponent')
    knee_freq, taus = convert_knee_val(fit_knee, fit_exp)

    # convert timescale into ms
    taus = taus * 1000

    # convert taus into a data frame
    taus_df = pd.DataFrame(taus.tolist(), index = colnames2use, columns = ["intsc"])
    taus_df = taus_df.T

    return(taus_df)


# boilerplate code to call main code for executing
if __name__ == '__main__':

    # parse arguments
    opts = parse_args()

    # main 4D time-series
    imgfile = opts.imgfile

    # output file
    outname = opts.outname

    # sampling rate
    srate = np.array(opts.srate, dtype = float)

    # load data
    data = load_data(imgfile)

    # compute intsc
    ntps = data.shape[0]
    res = compute_intsc(ts = data, numtps = ntps, srate = srate)

    # write result out to csv file
    outfile = "%s.csv" % outname
    export_csv = res.to_csv(outfile, index = None, header = True)

    print("Done")
