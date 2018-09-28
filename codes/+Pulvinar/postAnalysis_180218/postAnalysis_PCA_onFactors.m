%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/' ...
    'preAligned/multi-day_CoAoTdHoldRel_JanToMar/withGoodNeurons'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRel_JanToMar';

% add individual datasets
Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToMar/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180218', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180218;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
for nrun = 1:numel( par4 )
    rc2.addParams( par4( nrun ) );
end
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

return;

%% Post-running analysis - loading data and the output of LFADS
for nData = 1:length(dc.datasets)
    realData = dc.datasets(nData).loadData();
    r_real(nData) = R.Rstruct(realData.R);
end
% r_real = dc.datasets(1).loadData(); % get the original dataset (for all neurons)
% r_real = R.Rstruct(r_real.R); % put the dataset into R struct class




%% loading data and the output of LFADS
%for r_id = 1:length(rc2.runs)
factorsAllDays = [];
for r_id = 1
    run = rc2.runs(r_id); % pull out run information
    run.loadSequenceData(); % load sequence data in that run
    run.loadPosteriorMeans(); % load posterior mean in that run
    run.addPosteriorMeansToSeq();
    for nData = 1:length(run.sequenceData)
        currentDay = R.Rstruct(run.sequenceData{nData});
        factorsAllDays(nData).factors = [currentDay.r.factors];
        clear currentDay
    end
%     r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct
end


%% concatenate the factors to one matrix
factorsMatrix = [factorsAllDays.factors];
%% actually run PCA
[coeff, score, latent] = pca( factorsMatrix' );

%% calculate total variance explained
% normalize the latents to get variance explained by each PC
latent_tot = latent / sum(latent);
% calculate the cumulative variance explained
cum_variance = cumsum( latent_tot );
%% plot cum_variance

f1 = figure;
dimension = 1:40;
scatter(dimension, cum_variance)

%% Run params 
% run_1:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.2;   c_kl_co_weight = 0.2
% run_2:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_3:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.8;   c_kl_co_weight = 0.8
% run_4:    c_l2_gen_scale = 10;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_5:    c_l2_gen_scale = 50;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_6:    c_l2_gen_scale = 100;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5



%% Select the run and day you want to analyse
r_realCopy = r_real(6).copy();
r_lfadsWhole = RunID(1).r_lfads(6).copy();

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfadsWhole.r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfadsWhole.r(1).factors, 1);