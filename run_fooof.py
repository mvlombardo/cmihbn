# run_fooof.py

"""
Run FOOOF on EEG data. Will run FOOOF for each electrode.

Example usage:
python run_fooof.py --data RestingState_processed.txt --srate 500 --minfreq 3 --maxfreq 40
"""

# import modules
from fooof import FOOOF
import scipy as sp
import numpy as np
# import matplotlib.pyplot as plt
from optparse import OptionParser
import pandas as pd
from neurodsp import spectral


# function to parse input arguments
def parse_args():
    """
    Parse arguments.
    """
    parser=OptionParser()
    parser.add_option('--data',"",dest='data',help="RestingState_processed.txt file ex: --data RestingState_processed.txt",default=None)
    parser.add_option('--srate',"",dest='srate',help="Sampling rate file ex: --srate 500",default=None)
    parser.add_option('--minfreq',"",dest='minfreq',help="Minimum frequency file ex: --minfreq 3",default=3)
    parser.add_option('--maxfreq',"",dest='maxfreq',help="Maximum frequency file ex: --maxfreq 40",default=40)
    (options,args) = parser.parse_args()
    return(options)


# read in txt file with processed data as a data frame
def load_data(fname):
    data = pd.read_csv(fname, delimiter="\t")
    return(data)

# Calculate the psd
def calc_psd(data, srate):
    from scipy import signal
    (f, psd) = sp.signal.welch(data, srate, nperseg=1000)
    return(f, psd)

# Calculate psd using neurodsp.spectral.psd
def calc_spectral_psd(data, srate, method_type = "median"):
    (f, psd) = spectral.psd(x = data, Fs = srate, method = method_type, nperseg = srate*2)
    return(f, psd)


# boilerplate code to call main code for executing
if __name__ == '__main__':

    # Parse arguments
    opts = parse_args()
    datafile = opts.data
    srate = np.array(opts.srate, dtype = int)
    minfreq = np.array(opts.minfreq, dtype = int)
    maxfreq = np.array(opts.maxfreq, dtype = int)

    # load data
    data = load_data(datafile)

    # get names of electrodes
    colnames = list(data)
    # omit first row of rownames and last row of nothing
    electrodes = colnames[1:-1]

    # loop over electrodes
    for electrode in electrodes:
        # calculate PSD
        # (f, psd) = calc_psd(data[electrode], srate)
        (f, psd) = calc_spectral_psd(data[electrode], srate)

        # Initialize FOOOF object
        fm = FOOOF()

        # Define frequency range across which to model the spectrum
        freq_range = [minfreq, maxfreq]

        # Model the power spectrum with FOOOF, and print out a report
        fm.report(f, psd, freq_range)

        # save report and results to json file
        fname2save = "fooof_%s" % (electrode)
        fm.save_report(file_name = fname2save)
        fm.save(file_name = fname2save, save_results = True)
