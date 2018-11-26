# sliding_oof.py

"""
Run FOOOF over sliding windows.

Example usage:
python sliding_oof.py --subid NDARAA075AMK --task RestingState
"""


# import modules
import os
from fooof import FOOOF
import scipy as sp
import numpy as np
from optparse import OptionParser
import pandas as pd
from neurodsp import spectral
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser, OptionGroup


# function to parse input arguments
def parse_args():
    """
    Parse arguments.
    """

    # main options
    parser=OptionParser()
    parser.add_option('--subid',"",dest='subid',help="Subject ID ex: --subid NDARAA075AMK",default=None)
    parser.add_option('--task',"",dest='task',help="Task ex: --task RestingState",default=None)

    # additional options probably best left as the default
    extopts = OptionGroup(parser, "Additional options")
    extopts.add_option('',"--rootdir",dest='rootdir',help="Root directory where data is stored ex: --rootdir /Users/mvlombardo/Dropbox/HBN",default="/Users/mvlombardo/Dropbox/HBN")
    extopts.add_option('',"--srate",dest='srate',help="Sampling rate in Hz ex: --srate 500",default=500)
    extopts.add_option('',"--winsize",dest='winsize',help="Window size in number of samples ex: --winsize 2048",default=2048)
    extopts.add_option('',"--winstep",dest='winstep',help="Window step size in number of samples ex: --winstep 100",default=100)
    extopts.add_option('',"--minfreq",dest='minfreq',help="Minimum frequency for FOOOF ex: --minfreq 1",default=1)
    extopts.add_option('',"--maxfreq",dest='maxfreq',help="Maximum frequency for FOOOF ex: --maxfreq 50",default=50)
    extopts.add_option('',"--minwidth",dest='minwidth',help="Minimum frequency width limit for FOOOF ex: --minwidth 1",default=1)
    extopts.add_option('',"--maxwidth",dest='maxwidth',help="Maximum frequency width limit for FOOOF ex: --maxwidth 12",default=12)
    extopts.add_option('',"--verbose",action="store_true",dest='verbose',help="Set verbosity ex: --verbose",default=False)
    parser.add_option_group(extopts)

    (options,args) = parser.parse_args()
    return(options)


# read in txt file with processed data as a data frame
def load_data(fname):
    data = pd.read_csv(fname, delimiter="\t")
    return(data)


# Calculate psd using neurodsp.spectral.compute_spectrum
def calc_spectral_psd(data, srate, method_type = "median"):
    (f, psd) = spectral.compute_spectrum(sig = data, fs = srate, method = method_type, nperseg = srate*2)
    return(f, psd)


# grab windowed_data
def get_window_data(data, start_index, window_size_samples):
    window_start = start_index
    window_end = window_start + window_size_samples
    windowed_data = data.iloc[window_start:window_end,:]
    return(windowed_data)


# run fooof
def fooof_me(f, psd, minfreq, maxfreq, minwidth, maxwidth):
    freq_range = [minfreq, maxfreq]
    fm = FOOOF(peak_width_limits=[minwidth, maxwidth])
    fm.fit(f, psd, freq_range)
    res = fm.get_results()
    return(res)


# make heatmap of correlation matrix over electrodes sliding 1/f time-series
def make_heatmap(data, out_name = None, colormin=0, color_max=1, cmap2use="bwr", dpi2use=300):
    df4plot = data.astype(float)
    hm_plot = sns.heatmap(df4plot.corr(), vmin=color_min, vmax=color_max,cmap=cmap2use)
    if out_name is not None:
        fig2save = hm_plot.get_figure()
        fig2save.savefig(out_name, dpi=dpi2use)


# boilerplate code to call main code for executing
if __name__ == '__main__':

    # Parse arguments
    opts = parse_args()
    subid = opts.subid
    task = opts.task
    VERBOSE = opts.verbose

    # other parameters that stay fixed for HBN data
    srate = int(opts.srate)
    window_size_samples = int(opts.winsize)
    window_step_samples = int(opts.winstep)
    window_size_seconds = window_size_samples/srate

    # FOOOF parameters
    minwidth = int(opts.minwidth)
    maxwidth = int(opts.maxwidth)
    minfreq = int(opts.minfreq)
    maxfreq = int(opts.maxfreq)

    # paths
    rootdir = opts.rootdir
    ppdir = "%s/preproc_data" % rootdir
    subpath = "%s/%s" % (ppdir, subid)
    taskdir = "%s/%s" % (subpath, task)
    outpath = "%s/sliding_oof/%s/%s" % (rootdir, subid, task)
    os.system("mkdir -p %s" % outpath)
    fname2use = "%s/%s_processed.txt" % (taskdir, task)

    # load in preprocessed data
    ppdata = load_data(fname2use)
    ppdata = ppdata.drop(labels = [' ','Unnamed: 19'], axis = 1)
    electrodes = ppdata.columns

    # make window
    ppdata_size = ppdata.shape
    last_window_start = ppdata_size[0] - window_size_samples
    start_idx = np.arange(start=0,step=window_step_samples,stop=last_window_start, dtype=int)

    # pre-allocate df_oof data frame to store results in
    df_oof = pd.DataFrame(index = start_idx, columns = electrodes)

    # loop over windows
    for window_number, start_index in enumerate(start_idx):

        windowed_data = get_window_data(ppdata, start_index, window_size_samples)

        if VERBOSE:
            print("working on window %d" % window_number)
            print("window start %d, window end %d" % (window_start, window_end))

        # loop over electrodes
        for electrode in electrodes:

            # calculate power spectral density
            (f, psd) = calc_spectral_psd(windowed_data[electrode], srate)

            # Model the power spectrum with FOOOF
            res = fooof_me(f=f, psd=psd, minfreq=minfreq, maxfreq=maxfreq, minwidth=minwidth, maxwidth=maxwidth)
            df_oof.loc[start_index, electrode] = res.background_params[1]

    # make a heatmap of the correlation matrix of electrodes sliding 1/f slope time-series
    plotfname2save = "%s/%s_sliding_oof_heatmap.png" % (outpath, subid)
    make_heatmap(data = df_oof, out_name = plotfname2save)

    # write out a csv file
    fname2save = "%s/%s_sliding_oof_results.csv" % (outpath, subid)
    df_oof.to_csv(fname2save)
