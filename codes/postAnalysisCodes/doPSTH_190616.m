%%
loadpath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
            'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/180614data_rm_highCorr/170311_cueOnArrayOnTargetDim_HoldRel.mat'];

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
day_id = 6;
run = rc2.runs(r_id);
r_lfads = olapChopped.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms, 'rates');
assembled_lfads = R.Rstruct(r_lfads);

%%
tot_spikes = out.counts;
% nPieces = size(tot_spikes, 1);
%spikes = zeros(size(trueSpikes, 1), nSpikes);
clear spikes;
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
%loadpath = ['/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/M20170127/Gratings/M20170127_PUL_9M-g3-g4-g5-evokedSpiking-v10.mat'];
%data = load(loadpath);
%UE = data.UE;
%clear data
UE = UEs{day_id};

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
arrayStart_hold = arrayStart(UE.isHoldTrial);

stopTimes = stopInds(UE.isHoldTrial) - startInds(UE.isHoldTrial);

%% get binned timing for LFADS rate
cueStart_lfads = round(cueStart/par.spikeBinMs);
arrayStart_lfads = round(arrayStart/par.spikeBinMs);
arrayStart_lfads_hold = round(arrayStart_hold/par.spikeBinMs);
dimStart_lfads = round(dimStart/par.spikeBinMs);


%% set up aligned trials for real data and lfads data

half_ms = 500;
olapChopped_hold = olapChopped.r.r(UE.isHoldTrial);
%olapChopped_hold = olapChopped_hold((stopTimes - dimStart) > 302);

% for real data
assembled_hold_real = R.Rstruct(assembled_real.r(UE.isHoldTrial));
%assembled_hold_real = R.Rstruct(assembled_hold_real.r((stopTimes - dimStart) > 302));
%dimStart_real = dimStart_real((stopTimes - dimStart_real) > 302);

% for lfads data
assembled_hold_lfads = R.Rstruct(assembled_lfads.r(UE.isHoldTrial));
%assembled_hold_lfads = R.Rstruct(assembled_hold_lfads.r((stopTimes - dimStart) > 302));

%dimStart_lfads = dimStart_lfads((stopTimes - dimStart) > 302);

%%% get arrayshape info
%loadpath = ['/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/M20170127/MUA_GRATINGS/M20170308_PUL_63M-g3-g4-g6-g7-g8-g9-evokedSpiking-v12.mat'];
%data = load(loadpath);
%UE2 = data.UE;
%clear data

%% Align the trials for real data

% for cueAlign and arrayAlign
for i = 1:length(assembled_real.r)
    assembled_real.r(i).cueAlign = assembled_real.r(i).spikes(:, (cueStart(i)-(300-1)):(cueStart(i)+700) );
    assembled_real.r(i).arrayAlign = assembled_real.r(i).spikes(:, (arrayStart(i)-(half_ms-1)):(arrayStart(i)+half_ms) );
end

for i = 1:length(assembled_hold_real.r)
    assembled_hold_real.r(i).dimAlign = assembled_hold_real.r(i).spikes(:, (dimStart_real(i)-(half_ms-1)):(dimStart_real(i)+half_ms) );
    assembled_hold_real.r(i).arrayAlign = assembled_hold_real.r(i).spikes(:, (arrayStart_hold(i)-(300-1)):(arrayStart_hold(i)+700) );
end


%% align the trials for LFADS rate
half_lfads_cueLeft = 300/par.spikeBinMs;
half_lfads_cueRight = 700/par.spikeBinMs;
half_lfads = half_ms/par.spikeBinMs;
for i = 1:length(assembled_lfads.r)
    assembled_lfads.r(i).cueAlign = assembled_lfads.r(i).rates(:, (cueStart_lfads(i)-(half_lfads_cueLeft-1)):(cueStart_lfads(i)+half_lfads_cueRight) );
    assembled_lfads.r(i).arrayAlign = assembled_lfads.r(i).rates(:, (arrayStart_lfads(i)-(half_lfads-1)):(arrayStart_lfads(i)+half_lfads) );
end

for i = 1:length(assembled_hold_lfads.r)
    assembled_hold_lfads.r(i).dimAlign = assembled_hold_lfads.r(i).rates(:, (dimStart_lfads(i)-(half_lfads-1)):(dimStart_lfads(i)+half_lfads) );
    assembled_hold_lfads.r(i).arrayAlign = assembled_hold_lfads.r(i).rates(:, (arrayStart_lfads_hold(i)-(half_lfads_cueLeft-1)):(arrayStart_lfads_hold(i)+half_lfads_cueRight) );
end

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(assembled_real.r); % get trial number
nTimesRaw = size(assembled_real.r(1).cueAlign, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(assembled_real.r(1).cueAlign, 1); % get neuron nubmer


%% compute avg neuron firing rates for different cue locations for real spiking

sigma = 5;
binsize = par.spikeBinMs;
[AvgFiringRate_InRF_cue_real,AvgFiringRate_OffRF_cue_real, SmoothedSpikingMatrix_cue_real, TrialIndexInRF_cue_real, TrialIndexOffRF_cue_real, nTimesLFADS] = computeFRForPSTH_1(assembled_real, 'cueAlign', sigma, binsize, olapChopped.r.r, 'real', 'smooth');
[AvgFiringRate_InRF_array_real,AvgFiringRate_OffRF_array_real, SmoothedSpikingMatrix_array_real, TrialIndexInRF_array_real, TrialIndexOffRF_array_real, nTimesLFADS] = computeFRForPSTH_1(assembled_hold_real, ...
                                                  'arrayAlign', sigma, binsize, olapChopped_hold, 'real', 'smooth');
[AvgFiringRate_InRF_dim_real,AvgFiringRate_OffRF_dim_real, SmoothedSpikingMatrix_dim_real, TrialIndexInRF_dim_real, TrialIndexOffRF_dim_real, nTimesLFADS] = computeFRForPSTH_1(assembled_hold_real, ...
                                                  'dimAlign', sigma, binsize, olapChopped_hold, 'real', 'smooth');

%% compute avg neuron firing rates for different cue locations for LFADS rate
[AvgFiringRate_InRF_cue_lfads,AvgFiringRate_OffRF_cue_lfads, SmoothedSpikingMatrix_cue_lfads, TrialIndexInRF_cue_lfads, TrialIndexOffRF_cue_lfads, nTimesLFADS] = computeFRForPSTH_1(assembled_lfads, ...
                                                  'cueAlign', sigma, binsize, olapChopped.r.r, 'lfads', 'raw');
[AvgFiringRate_InRF_array_lfads,AvgFiringRate_OffRF_array_lfads, SmoothedSpikingMatrix_array_lfads, TrialIndexInRF_array_lfads, TrialIndexOffRF_array_lfads, nTimesLFADS] = computeFRForPSTH_1(assembled_hold_lfads, ...
                                                  'arrayAlign', sigma, binsize, olapChopped_hold, 'lfads', 'raw');
[AvgFiringRate_InRF_dim_lfads,AvgFiringRate_OffRF_dim_lfads, SmoothedSpikingMatrix_dim_lfads, TrialIndexInRF_dim_lfads, TrialIndexOffRF_dim_lfads, nTimesLFADS] = computeFRForPSTH_1(assembled_hold_lfads, ...
                                                  'dimAlign', sigma, binsize, olapChopped_hold, 'lfads', 'raw');


%% select good neurons (v12) - new good neurons after May 21, 18: why new good neurons???
% still need to update the 170320 - 170407 *****

%nIndices = [1 4 6 7 8 11 12 13 14 16 17 19 20 21 23 24 25 26 27 29 30 31 32];
% % 170127
% nIndices = [1 2 4 5 6 7 8 9 12 14 15 16 17 19 20 23 24 25 26 29 30 31 32];
% nIndices = [1 2 4 5 6 7 8 9 12 14 15 16 17 19 20 24 25 26 29 30 31 32]; % high corr removed
% % 170130
% nIndices = [1 5 7 8 10 12 13 27 30 32];
% % 170201
% nIndices = [3 4 5 7 11 12 15 17 18 23 28 30];
% % 170211
%  nIndices = [9 15 18 20 24 28 30 31 32 34 37 38 39 42 43 44 45 46 47 48 52 54 55 58 60 62 63 64];
%nIndices = [9 15 18 20 24 28 30 31 32 34 37 39 43 44 45 46 47 48 52 54 55 58 60 62 64]; % high Corr removed
% % 170308
% nIndices = [2 3 5 8 9 11 13 17 20 22 25 26 32 33 36 38 42 43 50 54 55 56 57 58];
nIndices = [2 3 5 8 9 11 13 17 20 22 25 26 32 33 36 38 42 50 54 55 56 57]; % high corr removed
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
%% plots of avg neuron firing rate and avg LFADS rate for different cue locations， plus spiking rasters
%savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDa%y_CO_AO_TD_HoldRel_JanToApr/' ...
%    'postAnalysis/withExternalInput_20180614/PSTH/NoSmoothing/arrayOnset/170130/'];


%savedirOne = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/reRun180614_20190530/PSTH/arr%ayOnset/170308';
%if ~isdir( savedirOne )
%    mkdir( savedirOne );
%end

cd(savedirOne);
clear set
plottingPSTHs_1(AvgFiringRate_InRF_array_real,AvgFiringRate_OffRF_array_real, SmoothedSpikingMatrix_array_real, AvgFiringRate_InRF_array_lfads,AvgFiringRate_OffRF_array_lfads, SmoothedSpikingMatrix_array_lfads, ...
              TrialIndexInRF_array_real, TrialIndexOffRF_array_real, nTimesLFADS, 'arrayOnset', nIndices, savedirOne, 3/10, '-300ms', '+700ms'  )

%savedirTwo = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/reRun180614_20190530/PSTH/cue%Onset/170308';
%if ~isdir( savedirTwo )
%    mkdir( savedirTwo );
%end

cd(savedirTwo);
clear set
plottingPSTHs_1(AvgFiringRate_InRF_cue_real,AvgFiringRate_OffRF_cue_real, SmoothedSpikingMatrix_cue_real, AvgFiringRate_InRF_cue_lfads,AvgFiringRate_OffRF_cue_lfads, SmoothedSpikingMatrix_cue_lfads, ...
              TrialIndexInRF_cue_real, TrialIndexOffRF_cue_real, nTimesLFADS, 'cueOnset', nIndices, savedirTwo, 3/10, '-300ms', '+700ms'  )


%%
%run180614.AvgFiringRate_InRF_cue_real = AvgFiringRate_InRF_cue_real;
%run180614.AvgFiringRate_OffRF_cue_real = AvgFiringRate_OffRF_cue_real;
%run180614.SmoothedSpikingMatrix_cue_real = SmoothedSpikingMatrix_cue_real;
%run180614.AvgFiringRate_InRF_cue_lfads = AvgFiringRate_InRF_cue_lfads;
%run180614.AvgFiringRate_OffRF_cue_lfads = AvgFiringRate_OffRF_cue_lfads;
%run180614.SmoothedSpikingMatrix_cue_lfads = SmoothedSpikingMatrix_cue_lfads;
%run180614.TrialIndexInRF_cue_real = TrialIndexInRF_cue_real;
%run180614.TrialIndexOffRF_cue_real = TrialIndexOffRF_cue_real;
%run180614.nTimesLFADS = nTimesLFADS;
%run180614.nIndices = nIndices;


%saveDataDir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/PSTH/SFN/data/';
%cd(saveDataDir)
%saveName = 'run180614';
%save(saveName, 'run180614');