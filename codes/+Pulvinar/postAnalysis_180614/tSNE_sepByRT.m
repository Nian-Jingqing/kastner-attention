%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
    'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRelSepForAO_JanToMar';

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
rc2 = Pulvinar.RunCollection(runRoot, 'withExternalInput_20180614', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180614;

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
par.c_ext_input_dim = 6;
par.c_inject_ext_input_to_gen = true;
par.c_keep_ratio = 1;

% the following params make sure this is a pbt run
par.doPBT = true;
par.PBTscript = '/snel/home/fzhu23/Projects/Pulvinar/codes/+Pulvinar/drive_scripts/pbt_script_run_manager_180614.py';


rc2.addParams( par );
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));


%% change day and define Directories here
dayName = 170127;
dayStr = num2str(dayName);
day_id = 1;
r_id = 1;
realDataBasePath = '/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag'; % from the pre-processed spiking data
rawDataBasePath = '/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/';
realDataFilePattern = '_cueOnArrayOnTargetDim_HoldRel.mat';

%% load real data
realDataFileName = [dayStr, realDataFilePattern];
realDataFileDir = fullfile(realDataBasePath, realDataFileName);
olapChopped = load(realDataFileDir);
olapChopped = olapChopped.combinedData;

%% get real data
trueSpikes = [olapChopped.r.r.spikes];
nSpikes = size(trueSpikes, 2);
trial_time_ms = 500;
trial_olap_ms = 100;
out = olapChopped.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );

%% get LFADS data
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
subFolderName = ['M20', dayStr];
timingBaseDir = fullfile(rawDataBasePath, subFolderName, 'Gratings');
timingFilePattern = fullfile(timingBaseDir, '*.mat'); % Change to whatever pattern you need.
timingFiles = dir(timingFilePattern);
firstTimingFileName = timingFiles(1).name;
firstTimingFile = fullfile(timingBaseDir, firstTimingFileName);
data = load(firstTimingFile);
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

% and also pick up good neurons here

% ######################################################################

allDayIndices{1} = [1 4 6 7 8 11 12 13 14 16 17 19 20 21 23 24 25 26 27 29 30 31 32];
% 170127
allDayIndices{2} = [1 2 4 5 6 7 8 9 12 14 15 16 17 19 20 23 24 25 26 29 30 31 32];
% 170130
allDayIndices{3} = [1 5 7 8 10 12 13 27 30 32];
% 170201
allDayIndices{4} = [3 4 5 7 11 12 15 17 18 23 28 30];
% 170211
allDayIndices{5} = [9 15 18 20 24 28 30 31 32 34 37 38 39 42 43 44 45 46 47 48 52 54 55 58 60 62 63 64];
% 170308
allDayIndices{6} = [2 3 5 8 9 11 13 17 20 22 25 26 32 33 36 38 42 43 50 54 55 56 57 58];
% 170311
nIndices = allDayIndices{day_id};

%% get RT and CL
rtVector = [olapChopped.r.r.rt];
rtVector_hold = [olapChopped_hold.rt];
clVector = [olapChopped.r.r.cueLoc];
clVector_hold = [olapChopped_hold.cueLoc];

%% ########### MODIFY THIS ####################
% define slow and faster RT
rtTagSlow = rtVector_hold > median(rtVector_hold)+0.05;
% any rt > median rt + 0.05s in rtVector is considered to be slow, indicated by
% logical 1
rtTagFast = rtVector_hold < median(rtVector_hold)-0.05;

% define timeLags that you want to analyze a
timeLagAll = -240:20:200;
% define sigma for smoothing real spiking data
sigma = 50;
% define the savedir to store the plots
savedir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/tSNE/170127';
% #############################################

%% smooth the real spiking data for tSNE
assembled_hold_real.smoothFieldInR( 'dimAlign', 'spike_smoothed', sigma, 1);
rbinned = assembled_hold_real.binData({'spike_smoothed'}, [par.spikeBinMs]);

%% normalize the data
for i = 1:length(assembled_hold_lfads.r)
    assembled_hold_lfads.r(i).normDimAlign = normalize(assembled_hold_lfads.r(i).dimAlign, 'centered');
end

for i = 1:length(rbinned)
    rbinned(i).normSpike_smoothed = normalize(rbinned(i).spike_smoothed, 'centered');
end



%% Actually performing t-SNE
timeAlign = nTimesRaw/2;
nTimesLFADS = nTimesRaw/par.spikeBinMs;
timeAlign_LFADS = nTimesLFADS/2;
timeLagAll_LFADS = timeLagAll/par.spikeBinMs;
for nLag = 1:length(timeLagAll)
%for nLag = 1

    timeLag_LFADS = timeLagAll_LFADS(nLag);
    timeLag = timeLagAll(nLag);
    f1 = figure
    tSNE_computePlot(timeLag_LFADS, timeAlign_LFADS, rtVector_hold, rbinned, 'real', rtTagSlow, rtTagFast, savedir, par.spikeBinMs, f1)
    tSNE_computePlot(timeLag_LFADS, timeAlign_LFADS, rtVector_hold, assembled_hold_lfads.r, 'lfads', rtTagSlow, rtTagFast, savedir, par.spikeBinMs, f1)
    suptitle(['TargetDim ' num2str(timeLag, '%+0.0f') 'ms']);
    set(f1, 'Position', [166 317 1287 510]);
    cd(savedir)
    print(f1,['TargetDim ' num2str(timeLag, '%+0.0f') 'ms'], '-dpng');
    close
end