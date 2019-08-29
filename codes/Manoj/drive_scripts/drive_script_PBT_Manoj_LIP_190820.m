%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
%% Locate and specify the datasets
datasetPath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/data_masked_highCorr_rm/';
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multi-Days';

% add individual datasets
% add individual datasets
Pulvinar.Dataset(dc, '0218_v2.mat');
Pulvinar.Dataset(dc, '0306_v2.mat');
Pulvinar.Dataset(dc, '0311_v2.mat');
Pulvinar.Dataset(dc, '0314_v2.mat');
Pulvinar.Dataset(dc, '0406_v2.mat');
Pulvinar.Dataset(dc, '0425_v2.mat');
Pulvinar.Dataset(dc, '0502_v2.mat');


% Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'tc_PBT_190820', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20190820;

% this script defines the run params
% Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
%for nrun = 1:numel( par4 )
%    rc2.addParams( par4( nrun ) );
%end

par = Pulvinar.RunParams;
par.spikeBinMs = 20;
par.c_batch_size = 2; % must be < 1/5 of the min trial count
par.trainToTestRatio = 10;
par.useAlignmentMatrix = true; % use alignment matrices

% the following params are specified in pbt script run manager
par.c_factors_dim = 40; 
%par.c_in_factors_dim = 30;
% par.c_co_dim = 0; % number of units in controller
% par.c_gen_dim = 64; % number of units in generator RNN
% par.c_ic_enc_dim = 64; % number of units in encoder RNN
% par.c_ci_enc_dim = 64; % number of units in encoder for controller input
% par.c_con_dim = 64; % controller dimensionality
% par.c_l2_gen_scale = 10;
% par.c_l2_con_scale = 10;
% par.c_kl_ic_weight = 0.2;
% par.c_kl_co_weight = 0.2;
par.c_learning_rate_stop = 5e-4;
par.c_learning_rate_decay_factor = 0.95;
%par.c_kl_increase_steps = 1;
%par.c_l2_increase_steps = 1;
%par.c_ext_input_dim = 15;
%par.c_inject_ext_input_to_gen = true;
par.c_keep_ratio = 0.7;

% the following params make sure this is a pbt run
par.doPBT = true;
par.PBTscript = '/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/Manoj/drive_scripts/pbt_script_run_manager_190820.py';


rc2.addParams( par );
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

%% Prepare LFADS input

rc2.prepareForLFADS();

rc2.writePBTShellScript()

