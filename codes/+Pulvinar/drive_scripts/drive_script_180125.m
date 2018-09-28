%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
datasetPath = ['/snel/share/share/derived/kastner/data_processed/singleSession/M20170608_PUL_all-g2-g3-g4-evokedSpiking/' ...
    'preAligned/CueOnArrayOnTargetDim_HoldRel/datasets'];

%% Locate and specify the datasets
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'preAligned_CO_AO_TD_HoldRel_20170608';

% add individual datasets
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/runs'];
rc = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180125', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc.version = 20180125;

Pulvinar.definePulvinarRunParams;

return;

%% Prepare LFADS input

% generate all of the data files LFADS needs to run everything
rc.prepareForLFADS(1);

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
%rc.writeShellScriptRunQueue('display', 50, 'maxTasksSimultaneously', 4, 'gpuList', [0 1], 'virtualenv', 'tensorflow');
rc.writeShellScriptRunQueue('display', 9, 'maxTasksSimultaneously', 4, 'gpuList', [0 1]);

%% check LFADS runs
r1 = rc.runs( 1 );
log = LFADS.Interface.read_fitlog( fullfile( r1.pathLFADSOutput, 'fitlog.csv' ) );
recTrain = cellfun( @(x) str2num(x), log(:, 9));
recValid = cellfun( @(x) str2num(x), log(:, 10));