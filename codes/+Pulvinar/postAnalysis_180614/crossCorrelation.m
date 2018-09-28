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
dayName = 170130;
dayStr = num2str(dayName);
day_id = 2;
r_id = 1;
lfpBasePath = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/continuousOverlapChop/';
realDataBasePath = '/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag'; % from the pre-processed spiking data
rawDataBasePath = '/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/';
lfpFilePattern = '_cueOnArrayOnTargetDim_HoldRel_lfp.mat';
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

%% load lfp data
       
lfpDayBasePath = fullfile(lfpBasePath, dayStr);
lfpFileName = [dayStr, lfpFilePattern];
lfpDir = fullfile(lfpDayBasePath, lfpFileName);
lfp_data = load(lfpDir);

%% get rid of outlier trials
lfp_noOutlier = lfp_data.r.r(~[lfp_data.r.r.isTrialOutlier]);

%% get rid of the common average channel
for i = 1:length(lfp_noOutlier)
    lfp_noOutlier(i).lfps = lfp_noOutlier(i).lfps(1:end - 1, :);
end

%% check if there are NaN in the LFP
for i = 1:length(lfp_noOutlier)
    a = isnan(lfp_noOutlier(i).lfps);
    b = any(a(:));
    if b
        disp('has NaN')
    end
end
        
    

%% bandpassFilter
filtHighCutoff = 20;
filtLowCutoff = 8;
Fs = 1000;

lfp_noOutlier = bandpassFilter_trialized( lfp_noOutlier, 'lfps', 'lfp_theta', filtHighCutoff, filtLowCutoff, Fs );

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

%% get cueOn, arrayOn, targetDim timing from real spiking data
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

timelag = 200;
cueStart = cueInds - startInds + timelag;
cueDelayEnd = arrayInds - startInds;
arrayStart = arrayInds - startInds + timelag;
dimStart = dimInds - startInds(UE.isHoldTrial);
% dimStart_real = dimStart;

% stopTimes = stopInds(UE.isHoldTrial) - startInds(UE.isHoldTrial);

%% get binned timing for LFADS rate
cueStart_lfads = round(cueStart/par.spikeBinMs);
cueDelayEnd_lfads = round(cueDelayEnd/par.spikeBinMs);
arrayStart_lfads = round(arrayStart/par.spikeBinMs);
dimStart_lfads = round(dimStart/par.spikeBinMs);

%% get rid of the outlier trials for spiking data, lfads data and timings
isTrialNotOutlier = ~[lfp_data.r.r.isTrialOutlier];
spiking_noOutlier = olapChopped.r.r(isTrialNotOutlier);
lfads_noOutlier = assembled_lfads.r(isTrialNotOutlier);
cueStart_noOutlier = cueStart(isTrialNotOutlier);
cueDelayEnd_noOutlier = cueDelayEnd(isTrialNotOutlier);
arrayStart_noOutlier = arrayStart(isTrialNotOutlier);
cueStart_lfads_noOutlier = cueStart_lfads(isTrialNotOutlier);
cueDelayEnd_lfads_noOutlier = cueDelayEnd_lfads(isTrialNotOutlier);
arrayStart_lfads_noOutlier = arrayStart_lfads(isTrialNotOutlier);
% find indices where are the no ouliers in the hold trials
indices_notOutlierInHold = isTrialNotOutlier(UE.isHoldTrial);
dimStart_noOutlier = dimStart(indices_notOutlierInHold);
dimStart_lfads_noOutlier = dimStart_lfads(indices_notOutlierInHold);
arrayStart_noOutlier_hold = arrayStart(isTrialNotOutlier' & UE.isHoldTrial);
arrayStart_lfads_noOutlier_hold = arrayStart_lfads(isTrialNotOutlier' & UE.isHoldTrial);

%{
%bothOutlierAndHoldTrial = false(size(cueStart)); % initialize bothOutlierAndHoldTrial
%bothOutlierAndHoldTrial([lfp_data.r.r.isTrialOutlier]' & UE.isHoldTrial) = true;
%makeOutlierInHoldToTwo = bothOutlierAndHoldTrial + UE.isHoldTrial;
%makeOutlierInHoldToTwo(makeOutlierInHoldToTwo == 0) = [];
%dimStart_noOutlier = dimStart(~(makeOutlierInHoldToTwo == 2));
%dimStart_lfads_noOutlier = dimStart_lfads(~(makeOutlierInHoldToTwo == 2));
%arrayStart_noOutlier_hold = arrayStart(~[lfp_data.r.r.isTrialOutlier]' & UE.isHoldTrial);
%arrayStart_lfads_noOutlier_hold = arrayStart_lfads(~[lfp_data.r.r.isTrialOutlier]' & UE.isHoldTrial);
%}

%% turn lfp_noOutlier to Rstruct and Resample LFP
lfp_noOutlier_rstruct = R.Rstruct(lfp_noOutlier);
rbinned_lfp_noOutlier(length(lfp_noOutlier_rstruct.r)).lfp_theta = 0;
for i = 1: length(lfp_noOutlier_rstruct.r)
    rbinned_lfp_noOutlier(i).lfp_theta = (resample(lfp_noOutlier_rstruct.r(i).lfp_theta', 1,4))';
end
% rbinned_lfp_noOutlier = lfp_noOutlier_rstruct.binData({'lfp_theta'}, [par.spikeBinMs]);

%% chop out first delay period for real spiking data and unbinned lfp_theta
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

% #######################################################################

unBinned_cueDelay(length(lfp_noOutlier)).lfp = 0; 
for i = 1:length(lfp_noOutlier)
    unBinned_cueDelay(i).lfp = lfp_noOutlier(i).lfp_theta(:, cueStart_noOutlier(i):cueDelayEnd_noOutlier(i) );
    unBinned_cueDelay(i).lfp = unBinned_cueDelay(i).lfp(nIndices, :);
    unBinned_cueDelay(i).spiking = spiking_noOutlier(i).spikes(:, cueStart_noOutlier(i):cueDelayEnd_noOutlier(i) );
    unBinned_cueDelay(i).spiking = full(unBinned_cueDelay(i).spiking);
end




%% chop out first delay period for lfads rates and binned lfp_theta
binned_cueDelay(length(lfp_noOutlier)).lfp = 0; 
for i = 1:length(lfp_noOutlier)
    binned_cueDelay(i).lfp = rbinned_lfp_noOutlier(i).lfp_theta(:, cueStart_lfads_noOutlier(i):cueDelayEnd_lfads_noOutlier(i) );
    binned_cueDelay(i).lfp = binned_cueDelay(i).lfp(nIndices, :);
    binned_cueDelay(i).lfadsRate = lfads_noOutlier(i).rates(:, cueStart_lfads_noOutlier(i):cueDelayEnd_lfads_noOutlier(i) );
end

%% chop out second delay period for real spiking data and unbinned lfp_theta
% first, need to find out which noOutlier is hold trial
indices_holdInNoOutlier = UE.isHoldTrial(isTrialNotOutlier);
% then let get rid of the release trials in the noOutlier trials
lfp_noOutlier_hold = lfp_noOutlier(indices_holdInNoOutlier);
spiking_noOutlier_hold = spiking_noOutlier(indices_holdInNoOutlier);
rbinned_lfp_noOutlier_hold = rbinned_lfp_noOutlier(indices_holdInNoOutlier);
lfads_noOutlier_hold = lfads_noOutlier(indices_holdInNoOutlier);


unBinned_arrayDelay(length(lfp_noOutlier_hold)).lfp = 0; 
for i = 1:length(lfp_noOutlier_hold)
    unBinned_arrayDelay(i).lfp = lfp_noOutlier_hold(i).lfp_theta(:, arrayStart_noOutlier_hold(i):dimStart_noOutlier(i) );
    unBinned_arrayDelay(i).lfp = unBinned_arrayDelay(i).lfp(nIndices, :);
    unBinned_arrayDelay(i).spiking = spiking_noOutlier_hold(i).spikes(:, arrayStart_noOutlier_hold(i):dimStart_noOutlier(i) );
    unBinned_arrayDelay(i).spiking = full(unBinned_arrayDelay(i).spiking);
end



%% chop out second delay period for lfads rates and binned lfp_theta
binned_arrayDelay(length(lfp_noOutlier_hold)).lfp = 0; 
for i = 1:length(lfp_noOutlier_hold)
    binned_arrayDelay(i).lfp = rbinned_lfp_noOutlier_hold(i).lfp_theta(:, arrayStart_lfads_noOutlier_hold(i):dimStart_lfads_noOutlier(i) );
    binned_arrayDelay(i).lfp = binned_arrayDelay(i).lfp(nIndices, :);
    binned_arrayDelay(i).lfadsRate = lfads_noOutlier_hold(i).rates(:, arrayStart_lfads_noOutlier_hold(i):dimStart_lfads_noOutlier(i) );
end

%% tmp test
for i = 1:length(lfp_noOutlier_hold)
    unit.lfp(i,:) = binned_arrayDelay(i).lfp(19,1:111);
    unit.lfads(i,:) = binned_arrayDelay(i).lfadsRate(19,1:111);
    unit.spikes(i,:) = unBinned_arrayDelay(i).spiking(19, 1:111*4);
end
%% tmp plot
f1 = figure;
s1 = subplot(3,1,1);
imagesc(unit.spikes);
s2 = subplot(3, 1, 2);
imagesc(unit.lfp)
s3 = subplot(3,1,3);
imagesc(unit.lfads);

%% select an equal time period for all trials for cueDelay and arrayDelay
cueDelayDuration = min(abs(cueStart_noOutlier - cueDelayEnd_noOutlier));
cueDelayDuration_lfads = min(abs(cueStart_lfads_noOutlier - cueDelayEnd_lfads_noOutlier));
arrayDelayDuration = min(abs(arrayStart_noOutlier_hold - dimStart_noOutlier));
arrayDelayDuration_lfads = min(abs(arrayStart_lfads_noOutlier_hold - dimStart_lfads_noOutlier));



%%
nNeurons = size(unBinned_cueDelay(1).lfp, 1);
nTrials = length(unBinned_cueDelay);
nTrials_hold = length(unBinned_arrayDelay);
cueLoc_list = unique(UE.cueLoc);
nCueLoc = length(cueLoc_list);
cueLoc_noOutlier = UE.cueLoc(isTrialNotOutlier);
cueLoc_noOutlier_hold = UE.cueLoc(isTrialNotOutlier' & UE.isHoldTrial);
rfLoc = olapChopped.r.r(1).rfloc;

%% cross-correlation for cueDelay and arrayDelay for different condition types
maxLag = 200;
xcorr_struct_cueDelay_inRF = compute_xcorr(nNeurons, unBinned_cueDelay, binned_cueDelay, cueDelayDuration, cueDelayDuration_lfads, maxLag, par, rfLoc, cueLoc_noOutlier, 'inRF');

xcorr_struct_cueDelay_oppositeRF = compute_xcorr(nNeurons, unBinned_cueDelay, binned_cueDelay, cueDelayDuration, cueDelayDuration_lfads, maxLag, par, rfLoc, cueLoc_noOutlier, 'oppositeRF');

xcorr_struct_arrayDelay_inRF = compute_xcorr(nNeurons, unBinned_arrayDelay, binned_arrayDelay, arrayDelayDuration, arrayDelayDuration_lfads, maxLag, par, rfLoc, cueLoc_noOutlier_hold, 'inRF');

xcorr_struct_arrayDelay_oppositeRF = compute_xcorr(nNeurons, unBinned_arrayDelay, binned_arrayDelay, arrayDelayDuration, arrayDelayDuration_lfads, maxLag, par, rfLoc, cueLoc_noOutlier_hold, 'oppositeRF');

%% plotting the cross-correlegram for cueDelay period

savedirBase = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/crossCorrelation/lfpFreq8To20_200msTimelag/testsepByCondition/cueDelay/'];
savedirOne = fullfile(savedirBase, dayStr);
clear set
for n = 1:nNeurons
    f1 = figure;
    sp1 = subplot(2,2,1);
    plottingXCorr(n, xcorr_struct_cueDelay_inRF, 'InRF', '-r', 'r_spikingCorrelation', 'Spikes', maxLag);
    hold on
    sp3 = subplot(2,2,3);
    plottingXCorr(n, xcorr_struct_cueDelay_inRF, 'InRF', '-b', 'r_lfadsCorrelation', 'LFADS Rates', maxLag);
    hold on
    sp2 = subplot(2,2,2);
    plottingXCorr(n, xcorr_struct_cueDelay_oppositeRF, 'OppositeRF', '-r', 'r_spikingCorrelation', 'Spikes', maxLag);
    hold on
    sp4 = subplot(2,2,4);
    plottingXCorr(n, xcorr_struct_cueDelay_oppositeRF, 'OppositeRF', '-b', 'r_lfadsCorrelation', 'LFADS Rates', maxLag);
    suptitle(['Cue Delay for Multi-unit ' int2str(nIndices(n))]);
    set(f1, 'Position', [229 79 1573 887]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');
    close
end

%% plotting the cross-correlegram for arrayDelay period

savedirBase = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/crossCorrelation/lfpFreq8To20_200msTimelag/sepByCondition/arrayDelay/'];
savedirOne = fullfile(savedirBase, dayStr);
clear set
for n = 1:nNeurons
    f1 = figure;
    sp1 = subplot(2,2,1);
    plottingXCorr(n, xcorr_struct_arrayDelay_inRF, 'InRF', '-r', 'r_spikingCorrelation', 'Spikes', maxLag);
    hold on
    sp3 = subplot(2,2,3);
    plottingXCorr(n, xcorr_struct_arrayDelay_inRF, 'InRF', '-b', 'r_lfadsCorrelation', 'LFADS Rates', maxLag);
    hold on
    sp2 = subplot(2,2,2);
    plottingXCorr(n, xcorr_struct_arrayDelay_oppositeRF, 'OppositeRF', '-r', 'r_spikingCorrelation', 'Spikes', maxLag);
    hold on
    sp4 = subplot(2,2,4);
    plottingXCorr(n, xcorr_struct_arrayDelay_oppositeRF, 'OppositeRF', '-b', 'r_lfadsCorrelation', 'LFADS Rates', maxLag);
    suptitle(['Array Delay for Multi-unit ' int2str(nIndices(n))]);
    set(f1, 'Position', [229 79 1573 887]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');
    close
end

%%    
    shadedErrorBar([], xcorr_struct_cueDelay_inRF(n).r_spikingCorrelation, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay_inRF(n).r_spikingCorrelation, 1)) }, 'lineProps', '-r');
    shadedErrorBar([], xcorr_struct_cueDelay_inRF(n).r_spikingCorrelation_shuffled, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay_inRF(n).r_spikingCorrelation_shuffled, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});
    set(gca,'XTick',[1 0.5*size(xcorr_struct_cueDelay(n).r_spikingCorrelation, 2) size(xcorr_struct_cueDelay(n).r_spikingCorrelation, 2)]);
    set(gca,'XTickLabels',{'-200','0','+200'});
    set(gca,'XLim',[0 size(xcorr_struct_cueDelay(n).r_spikingCorrelation, 2)])
    ylabel('cross-correlation');
    title('Cross-correlation Between Spikes and LFP')
    hold on
    
    sp2 = subplot(2,1,2);
    shadedErrorBar([], xcorr_struct_cueDelay(n).r_lfadsCorrelation, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 1)) }, 'lineProps', '-b');
    shadedErrorBar([], xcorr_struct_cueDelay(n).r_lfadsCorrelation_shuffled, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay(n).r_lfadsCorrelation_shuffled, 1)) }, 'lineProps', {'Color', [0.6, 0.6, 0.6]});
    set(gca,'XTick',[1 0.5*size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 2) size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 2)]);
    set(gca,'XTickLabels',{'-200','0','+200'});
    set(gca,'XLim', [0 size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 2)])
    ylabel('cross-correlation');
    xlabel('Time lag (ms)');
    title('Cross-correlation Between LFADS Rates and LFP')
    
    suptitle(['Cue Delay for Multi-unit ' int2str(nIndices(n))]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');
    close
end

%%
% ##################### no separation between cue locations #############################
%% cross-correlation for spikes vs LFP for the cue delay
maxLag = 200
allTrialIndices = 1: nTrials;
for n = 1:nNeurons
    for itrial = 1: nTrials
        spiking = unBinned_cueDelay(itrial).spiking(n,:);
        allOtherTrialIndices = allTrialIndices(allTrialIndices~=itrial);
        oneRandomTrial = randsample(allOtherTrialIndices,1);
        lfp_spiking = unBinned_cueDelay(itrial).lfp(n,:);
        lfp_spiking_shuffled = unBinned_cueDelay(oneRandomTrial).lfp(n,:);
        xcorr_struct_cueDelay(n).r_spikingCorrelation(itrial,:) = xcorr(spiking, lfp_spiking, maxLag);
        xcorr_struct_cueDelay(n).r_spikingCorrelation_shuffled(itrial,:) = xcorr(spiking, lfp_spiking_shuffled, maxLag);

        lfadsRate = binned_cueDelay(itrial).lfadsRate(n,:);
        lfp_lfads = binned_cueDelay(itrial).lfp(n,:);
        lfp_lfads_shuffled = binned_cueDelay(oneRandomTrial).lfp(n,:);
        xcorr_struct_cueDelay(n).r_lfadsCorrelation(itrial,:) = xcorr(lfadsRate, lfp_lfads, maxLag/par.spikeBinMs);
        xcorr_struct_cueDelay(n).r_lfadsCorrelation_shuffled(itrial,:) = xcorr(lfadsRate, lfp_lfads_shuffled, maxLag/par.spikeBinMs);
    end
end
    

%% cross-correlation for spikes vs LFP for the array delay
xcorr_struct_arrayDelay = struct;
% xcorr_struct(nNeurons).r_lfadsCorrelation = 0;
maxLag = 200;
allHoldTrialIndices = 1:nTrials_hold;
for n = 1:nNeurons
    for itrial = 1: nTrials_hold
        allOtherTrialIndices = allHoldTrialIndices(allHoldTrialIndices~=itrial);
        oneRandomTrial = randsample(allOtherTrialIndices,1);
        lfp_spiking = unBinned_arrayDelay(itrial).lfp(n,:);
        lfp_spiking_shuffled = unBinned_arrayDelay(oneRandomTrial).lfp(n,:);
        spiking = unBinned_arrayDelay(itrial).spiking(n,:);
        xcorr_struct_arrayDelay(n).r_spikingCorrelation(itrial,:) = xcorr(spiking, lfp_spiking, maxLag);
        xcorr_struct_arrayDelay(n).r_spikingCorrelation_shuffled(itrial,:) = xcorr(spiking, lfp_spiking_shuffled, maxLag);
        
        lfp_lfads = binned_arrayDelay(itrial).lfp(n,:);
        lfp_lfads_shuffled = binned_arrayDelay(oneRandomTrial).lfp(n,:);
        lfadsRate = binned_arrayDelay(itrial).lfadsRate(n,:);
        xcorr_struct_arrayDelay(n).r_lfadsCorrelation(itrial,:) = xcorr(lfadsRate, lfp_lfads, maxLag/par.spikeBinMs);
        xcorr_struct_arrayDelay(n).r_lfadsCorrelation_shuffled(itrial,:) = xcorr(lfadsRate, lfp_lfads_shuffled, maxLag/par.spikeBinMs);
    end
end

%% plotting the mean cross-correlation for each neuron
savedirBase = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/crossCorrelation/lfpFreq8To20_200msTimelag/cueDelay/'];
savedirOne = fullfile(savedirBase, dayStr);
clear set
for n = 1:nNeurons
    f1 = figure;
    sp1 = subplot(2,1,1);
    shadedErrorBar([], xcorr_struct_cueDelay(n).r_spikingCorrelation, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay(n).r_spikingCorrelation, 1)) }, 'lineProps', '-r');
    shadedErrorBar([], xcorr_struct_cueDelay(n).r_spikingCorrelation_shuffled, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay(n).r_spikingCorrelation_shuffled, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});
    set(gca,'XTick',[1 0.5*size(xcorr_struct_cueDelay(n).r_spikingCorrelation, 2) size(xcorr_struct_cueDelay(n).r_spikingCorrelation, 2)]);
    set(gca,'XTickLabels',{'-200','0','+200'});
    set(gca,'XLim',[0 size(xcorr_struct_cueDelay(n).r_spikingCorrelation, 2)])
    ylabel('cross-correlation');
    title('Cross-correlation Between Spikes and LFP')
    hold on
    
    sp2 = subplot(2,1,2);
    shadedErrorBar([], xcorr_struct_cueDelay(n).r_lfadsCorrelation, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 1)) }, 'lineProps', '-b');
    shadedErrorBar([], xcorr_struct_cueDelay(n).r_lfadsCorrelation_shuffled, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_cueDelay(n).r_lfadsCorrelation_shuffled, 1)) }, 'lineProps', {'Color', [0.6, 0.6, 0.6]});
    set(gca,'XTick',[1 0.5*size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 2) size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 2)]);
    set(gca,'XTickLabels',{'-200','0','+200'});
    set(gca,'XLim', [0 size(xcorr_struct_cueDelay(n).r_lfadsCorrelation, 2)])
    ylabel('cross-correlation');
    xlabel('Time lag (ms)');
    title('Cross-correlation Between LFADS Rates and LFP')
    
    suptitle(['Cue Delay for Multi-unit ' int2str(nIndices(n))]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');
    close
end

%% plotting the mean cross-correlation for each neuron for the array delay
savedirBase = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/crossCorrelation/lfpFreq8To20_200msTimelag/arrayDelay/'];
savedirOne = fullfile(savedirBase, dayStr);
clear set
for n = 1:nNeurons
    f1 = figure;
    sp1 = subplot(2,1,1);
    shadedErrorBar([], xcorr_struct_arrayDelay(n).r_spikingCorrelation, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_arrayDelay(n).r_spikingCorrelation, 1)) }, 'lineProps', '-r');
    shadedErrorBar([], xcorr_struct_arrayDelay(n).r_spikingCorrelation_shuffled, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_arrayDelay(n).r_spikingCorrelation_shuffled, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});
    set(gca,'XTick',[1 0.5*size(xcorr_struct_arrayDelay(n).r_spikingCorrelation, 2) size(xcorr_struct_arrayDelay(n).r_spikingCorrelation, 2)]);
    set(gca,'XTickLabels',{'-200','0','+200'});
    set(gca,'XLim',[0 size(xcorr_struct_arrayDelay(n).r_spikingCorrelation, 2)])
    ylabel('cross-correlation');
    title('Cross-correlation Between Spikes and LFP')
    hold on
    
    sp2 = subplot(2,1,2);
    shadedErrorBar([], xcorr_struct_arrayDelay(n).r_lfadsCorrelation, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_arrayDelay(n).r_lfadsCorrelation, 1)) }, 'lineProps', '-b');
    shadedErrorBar([], xcorr_struct_arrayDelay(n).r_lfadsCorrelation_shuffled, {@mean, @(x) std(x)./sqrt(size(xcorr_struct_arrayDelay(n).r_lfadsCorrelation_shuffled, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});
    set(gca,'XTick',[1 0.5*size(xcorr_struct_arrayDelay(n).r_lfadsCorrelation, 2) size(xcorr_struct_arrayDelay(n).r_lfadsCorrelation, 2)]);
    set(gca,'XTickLabels',{'-200','0','+200'});
    set(gca,'XLim', [0 size(xcorr_struct_arrayDelay(n).r_lfadsCorrelation, 2)])
    ylabel('cross-correlation');
    xlabel('Time lag (ms)');
    title('Cross-correlation Between LFADS Rates and LFP')
    
    suptitle(['Array Delay for Multi-unit ' int2str(nIndices(n))]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');
    close
end

