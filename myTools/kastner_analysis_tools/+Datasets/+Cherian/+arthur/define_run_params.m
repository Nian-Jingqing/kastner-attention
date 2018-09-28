%% Set parameters for the entire run collection


%% run 1: normal LFADS
par = cherian.RunParams;
par.doPBT = false;
%par.PBTscript = '/snel/home/mreza/projects/PBT_HP_opt/PBT_HP_opt/pbt_opt/pbt_script_run_manager.py';
%% Set parameters for the entire run collection

par = cherian.RunParams;
par.spikeBinMs = 4; % rebin the data at 2 ms
par.c_co_dim = 4; % no controller --> no inputs to generator
par.c_batch_size = 41; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; % and manually set it for multisession stitched models
par.useAlignmentMatrix = false; % use alignment matrices initial guess for multisession stitching

% *** MY DATASET PARAMETERS ***
par.mo_threshold = 0.25; % move onset threshold (percentage of max speed)
par.i_block = 1; % experiment block index from original data file
par.max_speed_min = 10; % minimum max speed of trial
par.pre_mo_ms = 200; % time before move onset
par.post_mo_ms = 700; % time after move onset
par.exclude_neuron_path = '/snel/home/lwimala/Projects/cherian-co/data-input/180402-arthur_robot053_s-neurons_to_remove.mat'; % path to data file containing 'all_neurons_to_remove' from x-corr analysis
par.dataset_type = 'center-out';                                                                                             
par.c_gen_dim = 128; % number of units in generator RNN
par.c_ic_enc_dim = 128; % number of units in encoder RNN

par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo

% add a single set of parameters to this run collection. Additional
% parameters can be added. LFADS.RunParams is a value class, unlike the other objects
% which are handle classes, so you can modify par freely.
rc.addParams(par);

%% run 2: normal LFADS w/ temporal spike jitter adjusted bin size
par = cherian.RunParams;
par.doPBT = false;
%par.PBTscript = '/snel/home/mreza/projects/PBT_HP_opt/PBT_HP_opt/pbt_opt/pbt_script_run_manager.py';
%% Set parameters for the entire run collection

par = cherian.RunParams;
par.spikeBinMs = 2; % rebin the data at 2 ms
par.c_co_dim = 4; % no controller --> no inputs to generator
par.c_batch_size = 41; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; % and manually set it for multisession stitched models
par.useAlignmentMatrix = false; % use alignment matrices initial guess for multisession stitching

% *** MY DATASET PARAMETERS ***
par.mo_threshold = 0.25; % move onset threshold (percentage of max speed)
par.i_block = 1; % experiment block index from original data file
par.max_speed_min = 10; % minimum max speed of trial
par.pre_mo_ms = 200; % time before move onset
par.post_mo_ms = 700; % time after move onse
par.dataset_type = 'center-out';                                                                                             
par.exclude_neuron_path = '/snel/home/lwimala/Projects/cherian-co/data-input/180402-arthur_robot053_s-neurons_to_remove.mat'; % path to data file containing 'all_neurons_to_remove' from x-corr analysis
                                                                                             
par.c_gen_dim = 128; % number of units in generator RNN
par.c_ic_enc_dim = 128; % number of units in encoder RNN

par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo


% altered parameters
par.c_temporal_spike_jitter_width = 2; % jittering spike times during training, in units of bin size
% add a single set of parameters to this run collection. Additional
% parameters can be added. LFADS.RunParams is a value class, unlike the other objects
% which are handle classes, so you can modify par freely.
rc.addParams(par);
