%% Set parameters for the entire run collection

%% run 1: normal LFADS w/ temporal spike jitter adjusted bin size

par = arthur.RunParams;
par.doPBT = false;
%par.PBTscript = '/snel/home/mreza/projects/PBT_HP_opt/PBT_HP_opt/pbt_opt/pbt_script_run_manager.py';
%% Set parameters for the entire run collection

par = arthur.RunParams;

% *** MY DATASET PARAMETERS ***

par.i_block = 1; % experiment block index from original data file
par.spikeBinMs = 2; % rebin the data at 2 ms
par.trialize_method = 'align-chop'; % method of splitting trials, 'align-chop' (typically center-out)  or 'overlap-chop' (random-walk)
par.dataset_type = 'center-out';
par.exclude_neuron_path = '/snel/home/lwimala/Projects/cherian-co/data-input/180402-arthur_robot053_s-neurons_to_remove.mat'; % path to data file containing 'all_neurons_to_remove' from x-corr analysis


% ***ALIGN-CHOP DATASET PARAMETERS

par.ap_threshold = 0.25; % threshold of max speed for alignment point of data
par.max_speed_min = 10; % minimum max speed of trial
par.pre_ap_ms = 200; % time (ms) before alignment point
par.post_ap_ms = 700; % time (ms) after alignment point

% ***OVERLAP-CHOP DATASET PARAMETERS

par.trial_time_ms = 500; % time (ms) of trial
par.trial_olap_ms = 100; % time (ms) of overlap between trials

% ***MODEL PARAMETERS

par.c_temporal_spike_jitter_width = 2; % jittering spike times during training, in units of bin size
par.c_gen_dim = 128; % number of units in generator RNN
par.c_ic_enc_dim = 128; % number of units in encoder RNN
par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo
par.c_co_dim = 4; % no controller --> no inputs to generator
par.c_batch_size = 41; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; % and manually set it for multisession stitched models
par.useAlignmentMatrix = false; % use alignment matrices initial guess for multisession stitching

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX RUN 1
rc.addParams(par);
