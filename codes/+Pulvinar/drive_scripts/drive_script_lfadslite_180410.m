%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
    'multi-unit/preAligned/multi-day_CoAoTdHoldRel_JanToApr/withGoodNeurons_issueFixed_chopped/'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRelSepForAO_JanToApr';

% add individual datasets
Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170320_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170324_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170327_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170329_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170331_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_lfadslite_20180410', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180410;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
par4( 1 ).c_in_factors_dim = 40;
par4( 1 ).c_keep_ratio = 0.5;
par4( 1 ).spikeBinMs = 4;
rc2.addParams( par4( 1 ) );
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

return;

%% Prepare LFADS input

% generate all of the data files LFADS needs to run everything
rc2.prepareForLFADS();

%%
% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
%rc.writeShellScriptRunQueue('display', 50, 'maxTasksSimultaneously', 4, 'gpuList', [0 1], 'virtualenv', 'tensorflow');
LFADSLITEPATH = '/snel/home/fzhu23/bin/lfadslite/run_lfadslite.py';
%rc2.writeShellScriptRunQueue('display', 9, 'maxTasksSimultaneously', 4, 'gpuList', [0 1], 'path_run_lfads_py', LFADSLITEPATH);
rc2.runs(1).writeShellScriptLFADSTrain('display', 9, 'cuda_visible_devices', [0], 'path_run_lfads_py', LFADSLITEPATH);