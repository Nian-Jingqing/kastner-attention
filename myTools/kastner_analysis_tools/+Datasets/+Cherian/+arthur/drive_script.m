%% Generate LFADS input datasets for PBT runs

addpath( '/snel/home/lwimala/Projects/lfads-run-manager/src' )
addpath( '/snel/home/lwimala/Projects/lfads/pbt-cherian' )
% build the dataset collection
datasetPath = '/snel/home/lwimala/Projects/lfads/pbt-cherian/data-input';

%% Locate and specify the datasets
dc = pbt_cherian.DatasetCollection(datasetPath);
dc.name = 'pbt-cherian-neural';

pbt_cherian.Dataset(dc, 'R-cherian-arthur-single-035.mat');

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = '/snel/home/lwimala/Projects/lfads/pbt-cherian/runs';
rc = pbt_cherian.RunCollection(runRoot, 'run002', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc.version = 20180302;

%% Set parameters for the entire run collection

par = pbt_cherian.RunParams;
par.spikeBinMs = 4; % rebin the data at 2 ms
par.c_co_dim = 0; % no controller --> no inputs to generator
par.c_batch_size = 256; % must be < 1/5 of the min trial count
par.c_factors_dim = 60; % and manually set it for multisession stitched models
par.useAlignmentMatrix = false; % use alignment matrices initial guess for multisession stitching
par.run_type = 'spikes';
par.data_switch = 0;
par.trial_time_ms = 500;
par.trial_olap_ms = 100;
par.c_gen_dim = 128; % number of units in generator RNN
par.c_ic_enc_dim = 128; % number of units in encoder RNN

par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo

% add a single set of parameters to this run collection. Additional
% parameters can be added. LFADS.RunParams is a value class, unlike the other objects
% which are handle classes, so you can modify par freely.
rc.addParams(par);

%% Add RunSpecs

% Run a single model for each dataset, and one stitched run with all datasets

% add each individual run
for iR = 1:dc.nDatasets
    runSpec = pbt_cherian.RunSpec(dc.datasets(iR).getSingleRunName(), dc, dc.datasets(iR).name);
    rc.addRunSpec(runSpec);
end


% adding a return here allows you to call this script to recreate all of
% the objects here for subsequent analysis after the actual LFADS models
% have been trained. The code below will setup the LFADS runs in the first
% place.

return;

%% Prepare LFADS input

% generate all of the data files LFADS needs to run everything
rc.prepareForLFADS();

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
rc.writeShellScriptRunQueue('display', 6, 'maxTasksSimultaneously', 4, 'gpuList', [0 1] );

%%

r = rc.runs( 1 );

r.loadSequenceData()
r.loadPosteriorMeans()
r.addPosteriorMeansToSeq()
