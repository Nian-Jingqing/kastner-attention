%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/fourLocations/'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRelSepForAO_JanToMar';

% add individual datasets
Pulvinar.Dataset(dc, '170529_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170601_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170608_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170613_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170615_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170618_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170621_cueOnArrayOnTargetDim_HoldRel.mat');
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
rc2 = Pulvinar.RunCollection(runRoot, 'withExternalInput_20180719', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180719;

% this script defines the run params
% Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
%for nrun = 1:numel( par4 )
%    rc2.addParams( par4( nrun ) );
%end

par = Pulvinar.RunParams;
par.spikeBinMs = 4;
par.c_batch_size = 2; % must be < 1/5 of the min trial count
par.trainToTestRatio = 50;
par.useAlignmentMatrix = true; % use alignment matrices

% the following params are specified in pbt script run manager
par.c_factors_dim = 30; 
par.c_in_factors_dim = 30;
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
par.c_kl_increase_steps = 1;
par.c_l2_increase_steps = 1;
par.c_ext_input_dim = 8;
par.c_inject_ext_input_to_gen = true;
par.c_keep_ratio = 1;

% the following params make sure this is a pbt run
par.doPBT = true;
par.PBTscript = '/snel/home/fzhu23/Projects/Pulvinar/codes/+Pulvinar/drive_scripts/pbt_script_run_manager_180719.py';


rc2.addParams( par );
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));


%% Post-running analysis - loading data and the output of LFADS
for nData = 1:length(dc.datasets)
    realData = dc.datasets(nData).loadData();
    r_real(nData) = R.Rstruct(realData.R);
end
% r_real = dc.datasets(1).loadData(); % get the original dataset (for all neurons)
% r_real = R.Rstruct(r_real.R); % put the dataset into R struct class




%% loading data and the output of LFADS
%for r_id = 1:length(rc2.runs)
for r_id = 1
    run = rc2.runs(r_id); % pull out run information
    run.loadSequenceData(); % load sequence data in that run
    run.loadPosteriorMeans(); % load posterior mean in that run
    run.addPosteriorMeansToSeq();
%     for nData = 1:length(run.sequenceData)
    nData = 5;
    RunID(r_id).r_lfads(nData) = R.Rstruct(run.sequenceData{nData});
%     end
%     r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct
end

%% Run params 
% run_1:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.2;   c_kl_co_weight = 0.2
% run_2:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_3:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.8;   c_kl_co_weight = 0.8
% run_4:    c_l2_gen_scale = 10;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_5:    c_l2_gen_scale = 50;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_6:    c_l2_gen_scale = 100;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5



%% Select the run and day you want to analyse
r_realCopy = r_real(5).copy();
r_lfadsWhole = RunID(1).r_lfads(5).copy();

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfadsWhole.r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfadsWhole.r(1).factors, 1);

%%
loadpath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/fourLocations/170613_cueOnArrayOnTargetDim_HoldRel.mat'];

olapChopped = load(loadpath);
olapChopped = olapChopped.combinedData;
%% get real data
trueSpikes = [olapChopped.r.r.spikes];
nSpikes = size(trueSpikes, 2);
trial_time_ms = 500;
trial_olap_ms = 100;
out = olapChopped.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );

%% get LFADS data
r_id = 1;
day_id = 4;
run = rc2.runs(r_id);
r_lfads = olapChopped.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms);
assembled_lfads = R.Rstruct(r_lfads);

%%
tot_spikes = out.counts;
% nPieces = size(tot_spikes, 1);
%spikes = zeros(size(trueSpikes, 1), nSpikes);
spikes(length(olapChopped.r.r)).spikes = 1;
%% assembly
nChops = 0;
for r_trial = 1: length(olapChopped.r.r)
    spikes(r_trial).spikes = zeros(size(olapChopped.r.r(r_trial).spikes));
    nPieces = ceil( (size( spikes(r_trial).spikes,2 ) - trial_time_ms ) / ( trial_time_ms - trial_olap_ms ) );
    for i = 1:nPieces
        spikesThisPiece = squeeze( tot_spikes( i+nChops, :, : ) );
        trialStart = 1 + (i-1)*trial_time_ms - (i-1)*trial_olap_ms;
        trialEnd = trialStart + trial_time_ms - 1;
%         if i~= size( tot_spikes, 1 )
%             trialStart = 1 + (i-1)*trial_time_ms - (i-1)*trial_olap_ms;
%             trialEnd = trialStart + trial_time_ms - 1;
%         else
%             trialStart = nSpikes - trial_time_ms + 1;
%             trialEnd = nSpikes;
%         end
        spikes(r_trial).spikes( :, trialStart:trialEnd ) = spikesThisPiece;
    end
    nChops = nChops + nPieces;
end

% %% trialize the continuous spikes
% r(length(olapChopped.r.r)).spikes = 1;
% past_nSpikes = 0;
% for i = 1: length(r)
%     nSpikesThisTrial = size(olapChopped.r.r(i).spikes,2);
%     start_ind = past_nSpikes +1;
%     end_ind = past_nSpikes + nSpikesThisTrial;
%     r(i).spikes = spikes(:, start_ind:end_ind);
%     past_nSpikes = end_ind;
% end
%%
assembled_real = R.Rstruct(spikes);

%% get timing info
loadpath = ['/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/fourLocations/M20170613/MUA_GRATINGS/M20170613_PUL_1M-g2-g3-g4-g5-evokedSpiking-v12.mat'];
data = load(loadpath);
UE = data.UE;
clear data

%% get cueOn, arrayOn, targetDim timing from real data
rawSampleRate = 1000;
sessStartTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1)*rawSampleRate;
% session start time - I set it to be the time when the monkey starts
% fixation in the first trial
sessEndTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
extraEndMs = 500;

startInds = round(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue * rawSampleRate - sessStartTime);
startInds(1) = 1;
stopInds = round(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice * rawSampleRate - sessStartTime);

cueInds = round(UE.cueOnset * rawSampleRate - sessStartTime);
arrayInds = round(UE.arrayOnset * rawSampleRate - sessStartTime);
dimInds = round(UE.targetDim * rawSampleRate - sessStartTime);

cueStart = cueInds - startInds;
arrayStart = arrayInds - startInds;
dimStart = dimInds - startInds(UE.isHoldTrial);
dimStart_real = dimStart;

stopTimes = stopInds(UE.isHoldTrial) - startInds(UE.isHoldTrial);

%% get binned timing for LFADS rate
cueStart_lfads = round(cueStart/par.spikeBinMs);
arrayStart_lfads = round(arrayStart/par.spikeBinMs);
dimStart_lfads = round(dimStart/par.spikeBinMs);


%% set up aligned trials for real data and lfads data

half_ms = 300;
olapChopped_hold = olapChopped.r.r(UE.isHoldTrial);
olapChopped_hold = olapChopped_hold((stopTimes - dimStart) > 302);

% for real data
assembled_hold_real = R.Rstruct(assembled_real.r(UE.isHoldTrial));
assembled_hold_real = R.Rstruct(assembled_hold_real.r((stopTimes - dimStart) > 302));
dimStart_real = dimStart_real((stopTimes - dimStart_real) > 302);

% for lfads data
assembled_hold_lfads = R.Rstruct(assembled_lfads.r(UE.isHoldTrial));
assembled_hold_lfads = R.Rstruct(assembled_hold_lfads.r((stopTimes - dimStart) > 302));

dimStart_lfads = dimStart_lfads((stopTimes - dimStart) > 302);


%% Align the trials for real data

% for cueAlign and arrayAlign
for i = 1:length(assembled_real.r)
    assembled_real.r(i).cueAlign = assembled_real.r(i).spikes(:, (cueStart(i)-(half_ms-1)):(cueStart(i)+half_ms) );
    assembled_real.r(i).arrayAlign = assembled_real.r(i).spikes(:, (arrayStart(i)-(half_ms-1)):(arrayStart(i)+half_ms) );
end

for i = 1:length(assembled_hold_real.r)
    assembled_hold_real.r(i).dimAlign = assembled_hold_real.r(i).spikes(:, (dimStart_real(i)-(half_ms-1)):(dimStart_real(i)+half_ms) );
end


%% align the trials for LFADS rate
half_lfads = 300/par.spikeBinMs;
for i = 1:length(assembled_lfads.r)
    assembled_lfads.r(i).cueAlign = assembled_lfads.r(i).rates(:, (cueStart_lfads(i)-(half_lfads-1)):(cueStart_lfads(i)+half_lfads) );
    assembled_lfads.r(i).arrayAlign = assembled_lfads.r(i).rates(:, (arrayStart_lfads(i)-(half_lfads-1)):(arrayStart_lfads(i)+half_lfads) );
end

for i = 1:length(assembled_hold_lfads.r)
    assembled_hold_lfads.r(i).dimAlign = assembled_hold_lfads.r(i).rates(:, (dimStart_lfads(i)-(half_lfads-1)):(dimStart_lfads(i)+half_lfads) );
end

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(assembled_real.r); % get trial number
nTimesRaw = size(assembled_real.r(1).cueAlign, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(assembled_real.r(1).cueAlign, 1); % get neuron nubmer


%% compute avg neuron firing rates for different cue locations for real spiking

sigma = 10;
binsize = par.spikeBinMs;
[AvgFiringRate_InRF_cue_real,AvgFiringRate_OffRF_cue_real, AvgFiringRate_PlusRF_cue_real, AvgFiringRate_MinusRF_cue_real, SmoothedSpikingMatrix_cue_real, TrialIndexInRF_cue_real, TrialIndexOffRF_cue_real, ...
 TrialIndexPlusRF_cue_real, TrialIndexMinusRF_cue_real, nTimesLFADS] = computeFRForPSTH_fourLocations(assembled_real, 'cueAlign', sigma, binsize, olapChopped.r.r, 'real');
[AvgFiringRate_InRF_array_real,AvgFiringRate_OffRF_array_real, AvgFiringRate_PlusRF_array_real, AvgFiringRate_MinusRF_array_real, SmoothedSpikingMatrix_array_real, TrialIndexInRF_array_real, ...
 TrialIndexOffRF_array_real, TrialIndexPlusRF_array_real, TrialIndexMinusRF_array_real, nTimesLFADS] = computeFRForPSTH_fourLocations(assembled_real, 'arrayAlign', sigma, binsize, olapChopped.r.r, 'real');
[AvgFiringRate_InRF_dim_real,AvgFiringRate_OffRF_dim_real, AvgFiringRate_PlusRF_dim_real, AvgFiringRate_MinusRF_dim_real, SmoothedSpikingMatrix_dim_real, TrialIndexInRF_dim_real, TrialIndexOffRF_dim_real, ...
 TrialIndexPlusRF_dim_real, TrialIndexMinusRF_dim_real, nTimesLFADS] = computeFRForPSTH_fourLocations(assembled_hold_real, 'dimAlign', sigma, binsize, olapChopped_hold, 'real');

%% compute avg neuron firing rates for different cue locations for LFADS rate
[AvgFiringRate_InRF_cue_lfads,AvgFiringRate_OffRF_cue_lfads, AvgFiringRate_PlusRF_cue_lfads, AvgFiringRate_MinusRF_cue_lfads, SmoothedSpikingMatrix_cue_lfads, TrialIndexInRF_cue_lfads, TrialIndexOffRF_cue_lfads, TrialIndexPlusRF_cue_lfads, TrialIndexMinusRF_cue_lfads, nTimesLFADS] = computeFRForPSTH_fourLocations(assembled_lfads, 'cueAlign', sigma, binsize, olapChopped.r.r, 'lfads');
[AvgFiringRate_InRF_array_lfads,AvgFiringRate_OffRF_array_lfads, AvgFiringRate_PlusRF_array_lfads, AvgFiringRate_MinusRF_array_lfads, SmoothedSpikingMatrix_array_lfads, TrialIndexInRF_array_lfads, ...
 TrialIndexOffRF_array_lfads, TrialIndexPlusRF_array_lfads, TrialIndexMinusRF_array_lfads, nTimesLFADS] = computeFRForPSTH_fourLocations(assembled_lfads, 'arrayAlign', sigma, binsize, olapChopped.r.r, 'lfads');
[AvgFiringRate_InRF_dim_lfads,AvgFiringRate_OffRF_dim_lfads, AvgFiringRate_PlusRF_dim_lfads, AvgFiringRate_MinusRF_dim_lfads, SmoothedSpikingMatrix_dim_lfads, TrialIndexInRF_dim_lfads, TrialIndexOffRF_dim_lfads, TrialIndexPlusRF_dim_lfads, TrialIndexMinusRF_dim_lfads, nTimesLFADS] = computeFRForPSTH_fourLocations(assembled_hold_lfads, 'dimAlign', sigma, binsize, olapChopped_hold, 'lfads');


%% select good neurons (v12) - new good neurons after May 21, 18: why new good neurons???
% still need to update the 170320 - 170407 *****

% nIndices = [1 4 6 7 8 11 12 13 14 16 17 19 20 21 23 24 25 26 27 29 30 31 32];
% % 170127
% nIndices = [1 2 4 5 6 7 8 9 12 14 15 16 17 19 20 23 24 25 26 29 30 31 32];
% % 170130
% nIndices = [1 5 7 8 10 12 13 27 30 32];
% % 170201
% nIndices = [3 4 5 7 11 12 15 17 18 23 28 30];
% % 170211
% nIndices = [9 15 18 20 24 28 30 31 32 34 37 38 39 42 43 44 45 46 47 48 52 54 55 58 60 62 63 64];
% % 170308
% nIndices = [2 3 5 8 9 11 13 17 20 22 25 26 32 33 36 38 42 43 50 54 55 56 57 58];
% % 170311
% nIndices = [4 5 9 12 16 20 24 25 29 30 31 34 35 37 38 40 41 44 49 50 51 52 53 54 56 57 59 61];
% % 170320
% nIndices = [1 4 7 8 11 12 13 14 15 16 18 21 25 27 28 29 30 33 34 36 39 42 43 44 45 47 48 49 50 51 52 57 58 60 61 62 63 64];
% % 170324
% nIndices = [1 2 3 4 5 10 11 13 15 16 18 19 21 22 24 25 27 31 32 35 39 40 41 42 43 44 45 47 50 52 53 54 55 56 57 58 59 60 61 62 63 64];
% % 170327
% nIndices = [1 2 3 5 6 7 8 12 14 16 17 18 20 22 25 26 27 28 29 31 32 34 38 39 40 41 47 48 49 50 51 52 53 54 55 56 60];
% % 170329 first half
% nIndices = [1 3 4 5 6 7 8 9 10 11 13 18 19 20 21 22 24 25 27 30 31 32 35 36 38 39 42 44 45 46 47 48 49 51 52 53 54 55 56 57 59 60 61 62 63 64];
% % 170329 second half
% nIndices = [7 9 10 11 13 14 17 24 25 29 30 31 32 34 36 37 41 46 47 50 51 52 55 57 58 59 62 64];
% % 170331 first half
% nIndices = [2 4 8 9 10 13 15 16 18 21 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 49 50 51 52 53 55 56 57 58 60 61 63 64];
% % 170331 second half
% nIndices = [16 19 20 23 24 32];
% % 170407
allDayIndices{1} = [7 8 9 14 17 18 19 24 29 31 35 36 49 50 52 53 56 60 62];
% 170529
allDayIndices{2} = [4 6 8 9 10 14 20 22 24 28 30 34 36 39 41 43 44 46 47 52 53 61 63];
% 170601
allDayIndices{3} = [1 2 3 4 17 22 28 31 33 36 37 43 46 49 55 59 63];
% 170608
allDayIndices{4} = [7 9 14 15 16 21 22 23 24 29 30 32 37 38 40 44 47 48 49 50 51 52 53 54 60];
% 170613
allDayIndices{5} = [1 3 5 7 8 10 12 13 16 20 24 30 36 40 42 45 47 48 53 54 59 63 64];
% 170615
allDayIndices{6} = [1 2 4 5 6 7 9 10 11 12 13 14 16 20 30 35 36 37 39 40 43 51 53 55];
% 170618
allDayIndices{7} = [2 4 7 8 12 14 16 20 22 38 39 43 44 46 48 54 56 62];
% 170621
nIndices = allDayIndices{day_id};

%% put all real and lfads average firing rate and spiking data into 1 struct for each aligned type

cueStruct_forPlotting.realRate_cond1 = AvgFiringRate_InRF_cue_real;
cueStruct_forPlotting.realRate_cond2 = AvgFiringRate_OffRF_cue_real;
cueStruct_forPlotting.realRate_cond3 = AvgFiringRate_PlusRF_cue_real;
cueStruct_forPlotting.realRate_cond4 = AvgFiringRate_MinusRF_cue_real;
cueStruct_forPlotting.realRaster = SmoothedSpikingMatrix_cue_real;
cueStruct_forPlotting.lfadsRate_cond1 = AvgFiringRate_InRF_cue_lfads;
cueStruct_forPlotting.lfadsRate_cond2 = AvgFiringRate_OffRF_cue_lfads;
cueStruct_forPlotting.lfadsRate_cond3 = AvgFiringRate_PlusRF_cue_lfads;
cueStruct_forPlotting.lfadsRate_cond4 = AvgFiringRate_MinusRF_cue_lfads;
cueStruct_forPlotting.lfadsRaster = SmoothedSpikingMatrix_cue_lfads;
cueStruct_forPlotting.trialIndexInRF = TrialIndexInRF_cue_real;
cueStruct_forPlotting.trialIndexOffRF = TrialIndexOffRF_cue_real;
cueStruct_forPlotting.trialIndexPlusRF = TrialIndexPlusRF_cue_real;
cueStruct_forPlotting.trialIndexMinusRF = TrialIndexMinusRF_cue_real;

arrayStruct_forPlotting.realRate_cond1 = AvgFiringRate_InRF_array_real;
arrayStruct_forPlotting.realRate_cond2 = AvgFiringRate_OffRF_array_real;
arrayStruct_forPlotting.realRate_cond3 = AvgFiringRate_PlusRF_array_real;
arrayStruct_forPlotting.realRate_cond4 = AvgFiringRate_MinusRF_array_real;
arrayStruct_forPlotting.realRaster = SmoothedSpikingMatrix_array_real;
arrayStruct_forPlotting.lfadsRate_cond1 = AvgFiringRate_InRF_array_lfads;
arrayStruct_forPlotting.lfadsRate_cond2 = AvgFiringRate_OffRF_array_lfads;
arrayStruct_forPlotting.lfadsRate_cond3 = AvgFiringRate_PlusRF_array_lfads;
arrayStruct_forPlotting.lfadsRate_cond4 = AvgFiringRate_MinusRF_array_lfads;
arrayStruct_forPlotting.lfadsRaster = SmoothedSpikingMatrix_array_lfads;
arrayStruct_forPlotting.trialIndexInRF = TrialIndexInRF_array_real;
arrayStruct_forPlotting.trialIndexOffRF = TrialIndexOffRF_array_real;
arrayStruct_forPlotting.trialIndexPlusRF = TrialIndexPlusRF_array_real;
arrayStruct_forPlotting.trialIndexMinusRF = TrialIndexMinusRF_array_real;

dimStruct_forPlotting.realRate_cond1 = AvgFiringRate_InRF_dim_real;
dimStruct_forPlotting.realRate_cond2 = AvgFiringRate_OffRF_dim_real;
dimStruct_forPlotting.realRate_cond3 = AvgFiringRate_PlusRF_dim_real;
dimStruct_forPlotting.realRate_cond4 = AvgFiringRate_MinusRF_dim_real;
dimStruct_forPlotting.realRaster = SmoothedSpikingMatrix_dim_real;
dimStruct_forPlotting.lfadsRate_cond1 = AvgFiringRate_InRF_dim_lfads;
dimStruct_forPlotting.lfadsRate_cond2 = AvgFiringRate_OffRF_dim_lfads;
dimStruct_forPlotting.lfadsRate_cond3 = AvgFiringRate_PlusRF_dim_lfads;
dimStruct_forPlotting.lfadsRate_cond4 = AvgFiringRate_MinusRF_dim_lfads;
dimStruct_forPlotting.lfadsRaster = SmoothedSpikingMatrix_dim_lfads;
dimStruct_forPlotting.trialIndexInRF = TrialIndexInRF_dim_real;
dimStruct_forPlotting.trialIndexOffRF = TrialIndexOffRF_dim_real;
dimStruct_forPlotting.trialIndexPlusRF = TrialIndexPlusRF_dim_real;
dimStruct_forPlotting.trialIndexMinusRF = TrialIndexMinusRF_dim_real;
%% plots of avg neuron firing rate and avg LFADS rate for different cue locations， plus spiking rasters
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180719/PSTH/NoSmoothing/arrayOnset/170613'];


%cd(savedirOne);
clear set
plottingPSTHs_fourLocations(arrayStruct_forPlotting, nTimesLFADS, 'arrayOnset', nIndices, savedirOne  )




    