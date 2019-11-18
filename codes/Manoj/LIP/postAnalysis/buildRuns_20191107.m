%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
%% Locate and specify the datasets
datasetPath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/newSignalProcessing/withPreAligned_withExtInp/lowerThresh/';
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = '10-days-goodChannels';

% add individual datasets
% add individual datasets
Pulvinar.Dataset(dc, '02082019_v2.mat');
Pulvinar.Dataset(dc, '02132019_v2.mat');
Pulvinar.Dataset(dc, '02142019_v2.mat');
Pulvinar.Dataset(dc, '02152019_v2.mat');
Pulvinar.Dataset(dc, '02182019_v2.mat');
Pulvinar.Dataset(dc, '02262019_v2.mat');
Pulvinar.Dataset(dc, '02282019_v2.mat');
Pulvinar.Dataset(dc, '03032019_v2.mat');
Pulvinar.Dataset(dc, '03062019_v2.mat');
Pulvinar.Dataset(dc, '03142019_v2.mat');


% Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'tc_newSP_PBT_191107_reRun', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 201911072;

% this script defines the run params
% Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
%for nrun = 1:numel( par4 )
%    rc2.addParams( par4( nrun ) );
%end

par = Pulvinar.RunParams;
par.spikeBinMs = 10;
par.c_batch_size = 2; % must be < 1/5 of the min trial count
par.trainToTestRatio = 10;
par.useAlignmentMatrix = true; % use alignment matrices

% the following params are specified in pbt script run manager
par.c_factors_dim = 40; 
%par.c_in_factors_dim = 40;
% par.c_co_dim = 4; % number of units in controller
                   %par.c_gen_dim = 64; % number of units in generator RNN
                    % par.c_ic_enc_dim = 64; % number of units in encoder RNN
% par.c_ci_enc_dim = 64; % number of units in encoder for controller input
% par.c_con_dim = 64; % controller dimensionality
% par.c_ic_dim = 64;
%par.c_l2_gen_scale = 1e-5;
%par.c_l2_con_scale = 1e-5;
%par.c_kl_ic_weight = 1e-5;
%par.c_kl_co_weight = 1e-5;
par.c_learning_rate_stop = 5e-4;
par.c_learning_rate_decay_factor = 0.95;
%par.c_kl_increase_epochs = 40;
%par.c_l2_increase_epochs = 40;
%par.c_n_epochs_early_stop = 50;
par.c_ext_input_dim = 15;
par.c_inject_ext_input_to_gen = true;
%par.c_keep_prob = 0.98;
par.c_keep_ratio = 0.7;

% the following params make sure this is a pbt run
par.doPBT = true;
par.PBTscript = '/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/Manoj/drive_scripts/pbt_script_run_manager_191107.py';


rc2.addParams( par );
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));