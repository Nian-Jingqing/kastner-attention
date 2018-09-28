%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/lorenz_example')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
addpath('/snel/home/fzhu23/bin/PBT_HP_opt/+PBT_analysis/make_pbt_run_plots.m')
%% Locate and specify the datasets
datasetPath = '/snel/home/fzhu23/Projects/lorenz_example/datasets';

% generate demo datasets
if ~exist(fullfile(datasetPath, 'dataset001.mat'), 'file')
    LFADS.Utils.generateDemoDatasets(datasetPath, 'nDatasets', 3);
end
dc = LorenzExperiment.DatasetCollection(datasetPath);
dc.name = 'lorenz_multiSession';

% add individual datasets
LorenzExperiment.Dataset(dc, 'dataset001.mat');
LorenzExperiment.Dataset(dc, 'dataset002.mat');
LorenzExperiment.Dataset(dc, 'dataset003.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/lorenz_runs/runs/test_multiday_PBT_180517';
rc2 = LorenzExperiment.RunCollection(runRoot, 'normalLFADS', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180517;

par = LorenzExperiment.RunParams;
par.name = 'multi-day'; % completely optional
par.useAlignmentMatrix = true; % use alignment matrices initial guess for multisession stitching

par.spikeBinMs = 5; % rebin the data at 5 ms
par.c_co_dim = 3; 
par.c_batch_size = 150; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; % and manually set it for multisession stitched models
par.c_in_factors_dim = 20;
par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_ci_enc_dim = 64; % number of units in encoder for controller input
par.c_con_dim = 64; % controller dimensionality
par.c_learning_rate_stop = 5e-4; 
par.c_learning_rate_decay_factor = 0.95;

% par.doPBT = true;
% par.PBTscript = '/snel/home/fzhu23/bin/PBT_HP_opt/pbt_opt/pbt_script_run_manager.py';



rc2.addParams( par );
rc2.addRunSpec(LorenzExperiment.RunSpec('all', dc, 1:dc.nDatasets));


%% Prepare LFADS input

% generate all of the data files LFADS needs to run everything
rc2.prepareForLFADS();

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
%rc.writeShellScriptRunQueue('display', 50, 'maxTasksSimultaneously', 4, 'gpuList', [0 1], 'virtualenv', 'tensorflow');
%rc2.writeShellScriptRunQueue('display', 9, 'maxTasksSimultaneously', 4, 'gpuList', [0 1]);
% rc2.writePBTShellScript()

LFADSPATH = '/snel/share/code/LFADS/research/lfads/run_lfads.py';
%rc2.writeShellScriptRunQueue('display', 9, 'maxTasksSimultaneously', 4, 'gpuList', [0 1], 'path_run_lfads_py', LFADSLITEPATH);
rc2.runs(1).writeShellScriptLFADSTrain('display', 9, 'cuda_visible_devices', [0], 'path_run_lfads_py', LFADSPATH);


%% analyze some PBT results

pbt_dir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/lorenz_runs/runs/test_multiday_PBT_180517/PBT/param_mFB3v1/all/pbt_run/';
results_save_dir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/lorenz_runs/postAnalysis/test_multiday_PBT_180517/PBT/PBT_plots/';
image_format = '-dpng';
make_pbt_run_plots( pbt_dir, results_save_dir, image_format )








