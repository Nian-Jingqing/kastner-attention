%% get real data
olapChopped = load(loadpath);
olapChopped = olapChopped.combinedData;
trueSpikes = [olapChopped.r.r.spikes];
nSpikes = size(trueSpikes, 2);
trial_time_ms = 500;
trial_olap_ms = 100;
out = olapChopped.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );

%% get LFADS data
run = rc2.runs(r_id);
r_lfads = olapChopped.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms, 'rates');
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
clear spikes
%% get timing info
rawDataBasePath = '/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/fourLocations/';
subFolderName = ['M20', dayStr];
timingBaseDir = fullfile(rawDataBasePath, subFolderName, 'MUA_GRATINGS');
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
%% savedir of plots of avg neuron firing rate and avg LFADS rate for different cue locations， plus spiking rasters
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20' runStr '/PSTH/NoSmoothing/' dayStr];

clear set
plottingPSTHs_fourLocations(arrayStruct_forPlotting, nTimesLFADS, 'arrayOnset', nIndices, savedirOne  )
plottingPSTHs_fourLocations(cueStruct_forPlotting, nTimesLFADS, 'cueOnset', nIndices, savedirOne  )
plottingPSTHs_fourLocations(dimStruct_forPlotting, nTimesLFADS, 'targetDim', nIndices, savedirOne  )
