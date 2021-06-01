function [EEG] = cmi_mat2set(cfg, subj_name, dataType)
% 
% CMI_MAT2SET convert CMI .mat files to .set
%   convert CMI .mat files into .set format for EEGlab
%

%% subject directories
subdir = fullfile(cfg.raw_data_dir,subj_name);
subdatadir = fullfile(subdir,dataType);
datafile2use = fullfile(subdatadir,sprintf('%s.mat',dataType));

%% load data
load(datafile2use, 'EEG')
eeg_raw = EEG;
sample_rate = eeg_raw.srate;

%% add channel info
chan_info_file = fullfile(cfg.code_dir,'GSN-HydroCel-129.sfp');
eeg_raw = pop_chanedit(eeg_raw, 'load', {chan_info_file,'filetype','sfp'});

%% add latency and latency_sec fields to eeg_raw.event(i_event)
n_event = length(eeg_raw.event);
event_cell = cell(n_event,1);

for i_event = 1:n_event
    
    event_cell{i_event,1} = eeg_raw.event(i_event).type;
    
    % latency
    eeg_raw.event(i_event).latency = eeg_raw.event(i_event).sample;
    
    % latency_sec
    eeg_raw.event(i_event).latency_sec = (eeg_raw.event(i_event).sample -1 )/ sample_rate;

end % for i_event = 1:n_event

%% make urevent field
eeg_raw = eeg_checkset(eeg_raw, 'makeur');

%% REMOVE junk/interval data
% BEFORE the first valid trigger
% and AFTER the last trigger
%    trigger1_sample = []; % ONSET
%    trigger2_sample = []; % OFFSET

if strcmp(dataType, 'RestingState') 
    % 90: start of Resting EEG paradigm
    % 20: eyes open start
    % 30: eyes closed start
    trigger_id = {'90  '; '20  '; '30  '};
elseif strcmp(dataType, 'DespicableMe')
    % 8_ = Start of Video: usually is video 3
    % 10_ = Stop of Video 
    trigger_id = [];
end % if strcmp(dataType, 'RestingState') 

sample_onset = []; 
sample_offset = [];

if strcmp(dataType, 'RestingState')
    
    % find first eyes open or closed epoch and use that as sample onset and
    % then find the last eyes open or closed epoch and use that as the
    % sample offset
    eyes_open_closed_event_idx = find(ismember(event_cell,trigger_id(2:3)));
    sample_onset = eeg_raw.event(eyes_open_closed_event_idx(1)).sample;
    sample_offset = eeg_raw.event(eyes_open_closed_event_idx(end)).sample;
    
    %elseif strcmp(file_name, 'desme.mat')
    % TO DO !!!
    
end % if strcmp(dataType, 'RestingState')

eeg_raw = pop_select(eeg_raw, 'point', [sample_onset-1,  sample_offset]);


%% save data to .set file
if cfg.do_save_eeg_set
    
    fprintf('... \n saving .set file')
    
    dir2save = fullfile(cfg.preproc_data_dir, subj_name, dataType);
    fname2save = fullfile(dir2save, sprintf('%s_%s.set',subj_name, dataType));
    unix(sprintf('mkdir -p %s',dir2save));
    pop_saveset(eeg_raw, 'filename', fname2save);
    
end % if do_save_eeg_set

% have eeg_raw be the new EEG structure
EEG = eeg_raw;

end % function cmi_mat2set

