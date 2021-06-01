% test_new_ppline.m

%% setup cfg structure
cfg.do_server = false;
cfg.project_dir = '/Users/mlombardo/Dropbox/cmi_test';
cfg.raw_data_dir = fullfile(cfg.project_dir,'data','raw');
cfg.preproc_data_dir = fullfile(cfg.project_dir,'data','preproc');
cfg.code_dir = fullfile(cfg.project_dir,'code');
cfg.eeglab_dir = '/Users/mlombardo/Dropbox/matlab/eeglab_20201226';
cfg.do_save_eeg_set = true;

%% install EEGlab
cd(cfg.eeglab_dir);
eeglab('nogui');
cd(cfg.code_dir);
addpath(genpath(cfg.code_dir));

%% set up preprocessing parameters
pp.downsample_rate = 250; % sampling rate to downsample to in Hz 

pp.hpf_cutoff = 1; % High pass filter in Hz 
pp.lpf_cutoff = 80; % low pass filter cutoff in Hz
%pp.lpf_cutoff = eeg_struct.srate/2 -1;

pp.line_noise_freq = 60; % line noise frequency in Hz 

pp.do_chan_rejection = true; % reject channels
% channels to reject
pp.chan_toreject = {
                'E127','E126',...
    'E25','E21','E17','E14','E8',...
    'E128','E32',          'E1','E125',...
    'E48','E43',            'E120','E119',...
    'E49',                          'E113',...
    'E56','E63',             'E99', 'E107',...
    'E68', 'E73', 'E81', 'E88', 'E94',...
               };
           
pp.do_chan_interp_pruning = true; % interpolate channels
% channels to interpolate
pp.chan_interp_prune =  {
    'E18','E10',...
    'E38','E121', ...
    'E44','E114', ...
    'E28','E117',...
    'E47','E98',...
    };

pp.do_waveletICA = false;    % <<< alternative a) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pp.do_cleanraw = true;      % cleanraw data + Artifact Subspace Reconstruction + avgref + ICA
    
pp.do_plot_chan = false;
pp.do_plot_PSD = true;
    
%pp.do_save_cleanline = false;
%pp.do_save_wavclean_ICA = true;
pp.do_save_cleanraw_avgref_ICA = true;

%pp.do_save_wavclean_nobadICA = true;
pp.do_save_cleanraw_avgref_nobadICA = false;

pp.do_save_fig = true;
pp.do_scorepoch = false;
pp.do_save_score = true;

% add pp to cfg
cfg.pp = pp;



%% subject ids
subids = {'NDARAM848GTE','NDARDH086ZKK'};

%% data types
dataTypes = {'RestingState'};

%% Full pipeline

% loop over subjects
for isub = 1:length(subids)
    
    sub_idx = isub;
    subid2use = subids{sub_idx};
    
    % loop over data types
    for idt = 1:length(dataTypes)
        
        dataType2use = dataTypes{idt};
        
        % Step 1: convert .mat to .set ==> cmi_mat2set.m
        [EEG] = cmi_mat2set(cfg, subid2use, dataType2use);
        cfg.setfilename{isub,1} = fullfile(cfg.preproc_data_dir, ...
            subid2use, ...
            dataType2use, ...
            sprintf('%s_%s.set', subid2use, dataType2use));
        
        % Step 2: Run preprocessing pipeline ==> cmi_preproc_land.m
        [EEG] = cmi_preproc_land(EEG, cfg, subid2use, sub_idx, dataType2use);
        
    end % for idt = 1:length(dataTypes)

end % for isub = 1:length(subids)
