% preproc_batch.m
% run preprocessing over subjects

% rootpath where underneath are individual subject directories
rootpath = '~/HBN';

% code path
codepath = fullfile(rootpath, 'code');

% channel loc file
channelLocFile = fullfile(codepath,'GSN-HydroCel-129_mvl_edit.sfp');

% sampling rate
samplingRate = 500;

% channels to use
channels2use = {'FP1','FP2','F3','F4','F7','F8','C3','C4','T3','T4', ...
    'PZ','O1','O2','T5','T6','P3','P4','Fz'};

% subjects to run
subids = {'NDARBM213BEA','NDARKA429EPF','NDARKX761BH9','NDARVH070WH6'};

% main loop over subjects
for isub = 1:length(subids)
    cd(codepath);
    datafile = fullfile(rootpath, subids{isub}, 'raw', 'RestingState.mat');
    preprocHAPPE(datafile, samplingRate, channels2use, channelLocFile);
end % for isub

cd(codepath);
