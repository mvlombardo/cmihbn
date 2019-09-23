"""
Compute intrinsic neural timescale (intsc) from time-series data.

Watanabe, T., Rees, G., & Masuda, N. (2019). Atypical intrinsic neural timescale in autism. eLife, 8, e42256.

Example usage

imgfile=~/data/rest_pp01.nii.gz
outfile=~/data/rest_pp01
srate=500
nlags=150000

python intsc_eeg.py -i $imgfile -o $outfile --srate 500 --nlags $nlags
"""

# import libraries
import numpy as np
import pandas as pd
import scipy as sp
import math
# import nibabel as nib
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
    parser.add_option('-v',"--verbose", \
                      dest='verbose', \
                      action='store_true', \
                      help="Use if you want verbose feedback while script runs.", \
                      default=False)
    parser.add_option('',"--srate", \
                      dest='srate',\
                      help="sampling rate in Hz", \
                      default=None)
    parser.add_option('',"--nlags", \
                      dest='nlags',\
                      help="number of lags, specified as number of samples", \
                      default=10000)

    (options,args) = parser.parse_args()

    return(options)


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


# function to find next power of 2
def nextpow2(x):
    """
    Get next power of 2
    """

    result = math.ceil(sp.log2(np.abs(x)))
    result = int(result)
    return(result)


# function to mean center data
def mean_center(data, n, verbose):
    """
    Mean center data
    """

    if verbose:
        print("Mean centering data")

    # compute the mean over time
    time_mean = np.nanmean(data)

    # take data and minus off the mean
    result = data - time_mean

    return(result)


# main function for computing intrinsic neural timescale
def compute_intsc(ts, numtps, nlags, srate):
    """
    Main function for computing intrinsic neural timescale
    """

    # figure out tr in seconds
    tr = 1/srate

    # find next power of 2
    p2use = nextpow2(numtps)

    # figure out n for FFT
    nFFT = 2**(p2use+1)

    # FFT data
    f = sp.fft(ts, n=nFFT)

    # multiply FFT result by its complex conjugate
    f = np.multiply(f, np.conjugate(f))

    # inverse FFT
    acf = sp.ifft(f)

    # grab nlags
    acf = acf[0:nlags+1]

    # normalize
    acf = np.divide(acf,acf[0])

    # convert to real numbers
    acf = np.real(acf)

    # sum positive acf
    positive_acf = np.nansum(acf[acf>0])-acf[0]

    # multiply by TR
    intsc = positive_acf * (srate/1000)

    return(intsc)


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

    # nlags
    nlags = np.array(opts.nlags, dtype = int)

    # verbose flag
    verbose = opts.verbose

    # load data
    data = load_data(imgfile)

    # loop over channels and compute intsc
    res = pd.DataFrame(index = range(1), columns = data.columns)
    for column in data.columns:
        # grab specific channel
        ts = data[column]

        # number of time points
        ntps = len(ts)

        # mean center time-series
        ts_mc = mean_center(data = ts, n = ntps, verbose = verbose)

        # compute_intsc
        res[column] = compute_intsc(ts = ts_mc, numtps = ntps, nlags = nlags, \
                                    srate = srate)

    # write result out to csv file
    outfile = "%s.csv" % outname
    export_csv = res.to_csv(outfile, index = None, header = True)

    print("Done")
