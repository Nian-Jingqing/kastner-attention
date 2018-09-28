%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
addpath(genpath('/snel/home/fzhu23/bin/chronux_2_12'))
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

%% get rid of the common average channel and drifting neurons
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
for i = 1:length(lfp_noOutlier)
    lfp_noOutlier(i).lfps = lfp_noOutlier(i).lfps(1:end - 1, :); % get rid of commond ave channel
    lfp_noOutlier(i).lfps = lfp_noOutlier(i).lfps(nIndices, :); % get rid of drifting neurons
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
filtHighCutoff = 80;
filtLowCutoff = 1;
Fs = 1000;
lfp_noOutlier = bandpassFilter_trialized( lfp_noOutlier, 'lfps', 'lfp_filtered', filtHighCutoff, filtLowCutoff, Fs );

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

cueStart = cueInds - startInds;
arrayStart_hold = arrayInds(UE.isHoldTrial) - startInds(UE.isHoldTrial);
dimStart = dimInds - startInds(UE.isHoldTrial);
% dimStart_real = dimStart;

% stopTimes = stopInds(UE.isHoldTrial) - startInds(UE.isHoldTrial);

%% get rid of the outlier trials for spiking data, lfads data and timings
isTrialNotOutlier = ~[lfp_data.r.r.isTrialOutlier];
spiking_noOutlier = olapChopped.r.r(isTrialNotOutlier);
cueStart_noOutlier = cueStart(isTrialNotOutlier);
% find indices where are the no ouliers in the hold trials
indices_notOutlierInHold = isTrialNotOutlier(UE.isHoldTrial);
dimStart_noOutlier = dimStart(indices_notOutlierInHold);
arrayStart_noOutlier_hold = arrayStart_hold(indices_notOutlierInHold);

%% retrieve rfLoc
rfLoc = olapChopped.r.r(1).rfloc;

%% set up parameters for preparing the neuronStruct for computing coherence 
slidingStart = -200;
slidingEnd = 300;
slidingStep = 50;
window_size = 300;
condType = 'OppositeRF';
alignType = 'arrayOnset';
if strcmp(alignType, 'cueOnset')
    spiking_struct = spiking_noOutlier;
    lfp_struct = lfp_noOutlier;
    alignInds = cueStart_noOutlier;
    cueLoc = UE.cueLoc(isTrialNotOutlier);
elseif strcmp(alignType, 'arrayOnset')
    indices_holdInNoOutlier = UE.isHoldTrial(isTrialNotOutlier);
    spiking_struct = spiking_noOutlier(indices_holdInNoOutlier);
    alignInds = arrayStart_noOutlier_hold;
    lfp_struct = lfp_noOutlier(indices_holdInNoOutlier);
    cueLoc = UE.cueLoc(isTrialNotOutlier' & UE.isHoldTrial);
end

%% Preparing the neuronStruct for computing coherence
[neuronStruct] = prepForWindowedCoherence(slidingStart, slidingEnd, slidingStep, window_size, spiking_struct, lfp_struct, alignInds, condType, cueLoc, rfLoc);

%% compute spike-field coherence with Chronux toolbox
% set up params
TW = 2; %time-bandwidth product
K = 3; % Number of tapers
params.tapers = [TW K];
params.pad = -1;
params.Fs = 1000;
params.fpass = [0 80];
params.trialave = 1;
nNeurons = length(neuronStruct);
for n = 1:length(neuronStruct)
    trialNum_thisNeuron = size(neuronStruct(n).window(1).spiking,1);
    for i = 1:length(neuronStruct(1).window)
        data1 = neuronStruct(n).window(i).lfp';
        data2 = neuronStruct(n).window(i).spiking';
        [C,phi,S12,S1,S2,f,zerosp] = coherencycpb(data1, data2, params);
        neuronStruct(n).window(i).C = atanh(C) - 1/(2*K*trialNum_thisNeuron - 2);
        neuronStruct(n).window(i).freq = f;
    end
end
%% assign neuronStruct to certain catagory
array_OppositeRF = neuronStruct;

%% plotting
savedirOne = '/snel/share/share/derived/kastner/nonLFADS_analysis/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/spike_field_coherence/170130/individualNeurons';
for n = 1: nNeurons
    f1 = figure;
    cue_InRF(n).coherence = flip([cue_InRF(n).window.C], 1);
    cue_OppositeRF(n).coherence = flip([cue_OppositeRF(n).window.C], 1);
    array_InRF(n).coherence = flip([array_InRF(n).window.C], 1);
    array_OppositeRF(n).coherence = flip([array_OppositeRF(n).window.C], 1);
    s1 = subplot(2,2,1)
    imagesc(cue_InRF(n).coherence)
    set(gca,'XTick',[1 4 11]);
    set(gca, 'XTickLabels', {'-150', '0', '+350'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s1, 'CueOnset - InRF');
    ylabel('Frequence (Hz)');

    s3 = subplot(2,2,3)
    imagesc(cue_OppositeRF(n).coherence)
    set(gca,'XTick',[1 4 11]);
    set(gca, 'XTickLabels', {'-150', '0', '+350'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s3, 'CueOnset - OppositeRF');
    ylabel('Frequence (Hz)');
    xlabel('Peri-cue time (ms)')

    s2 = subplot(2,2,2)
    imagesc(array_InRF(n).coherence)
    set(gca,'XTick',[1 5 10]);
    set(gca, 'XTickLabels', {'-200', '0', '+250'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s2, 'ArrayOnset - InRF');
    ylabel('Frequence (Hz)');

    s4 = subplot(2,2,4)
    imagesc(array_OppositeRF(n).coherence)
    set(gca,'XTick',[1 5 10]);
    set(gca, 'XTickLabels', {'-200', '0', '+250'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s4, 'ArrayOnset - OppositeRF');
    ylabel('Frequence (Hz)');
    xlabel('Peri-array time (ms)')

    suptitle(['Multi-unit ' int2str(nIndices(n))]);
    set(f1, 'Position', [452 136 1053 761]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');

    close;
end

%% compute all neuron avg
allNeuronAve_cue_InRF = zeros(size([cue_InRF(1).window.C]));
allNeuronAve_cue_OppositeRF = zeros(size([cue_OppositeRF(1).window.C]));
allNeuronAve_array_InRF = zeros(size([array_InRF(1).window.C]));
allNeuronAve_array_OppositeRF = zeros(size([array_OppositeRF(1).window.C]));
for n = 1: nNeurons
    cue_InRF(n).coherence = flip([cue_InRF(n).window.C], 1);
    cue_OppositeRF(n).coherence = flip([cue_OppositeRF(n).window.C], 1);
    array_InRF(n).coherence = flip([array_InRF(n).window.C], 1);
    array_OppositeRF(n).coherence = flip([array_OppositeRF(n).window.C], 1);
    
    allNeuronAve_cue_InRF = allNeuronAve_cue_InRF + cue_InRF(n).coherence;
    allNeuronAve_cue_OppositeRF = allNeuronAve_cue_OppositeRF + cue_OppositeRF(n).coherence;
    allNeuronAve_array_InRF = allNeuronAve_array_InRF + array_InRF(n).coherence;
    allNeuronAve_array_OppositeRF = allNeuronAve_array_OppositeRF + array_OppositeRF(n).coherence;
end
allNeuronAve_cue_InRF = allNeuronAve_cue_InRF/nNeurons;
allNeuronAve_cue_OppositeRF = allNeuronAve_cue_OppositeRF/nNeurons;
allNeuronAve_array_InRF = allNeuronAve_array_InRF/nNeurons;
allNeuronAve_array_OppositeRF = allNeuronAve_array_OppositeRF/nNeurons;

%%
savedirOne = '/snel/share/share/derived/kastner/nonLFADS_analysis/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/spike_field_coherence/170130/allNeuronAvg';
f1 = figure;
s1 = subplot(2,2,1)
    imagesc(allNeuronAve_cue_InRF)
    set(gca,'XTick',[1 4 11]);
    set(gca, 'XTickLabels', {'-150', '0', '+350'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s1, 'CueOnset - InRF');
    ylabel('Frequence (Hz)');

    s3 = subplot(2,2,3)
    imagesc(allNeuronAve_cue_OppositeRF)
    set(gca,'XTick',[1 4 11]);
    set(gca, 'XTickLabels', {'-150', '0', '+350'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s3, 'CueOnset - OppositeRF');
    ylabel('Frequence (Hz)');
    xlabel('Peri-cue time (ms)')

    s2 = subplot(2,2,2)
    imagesc(allNeuronAve_array_InRF)
    set(gca,'XTick',[1 5 10]);
    set(gca, 'XTickLabels', {'-200', '0', '+250'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s2, 'ArrayOnset - InRF');
    ylabel('Frequence (Hz)');

    s4 = subplot(2,2,4)
    imagesc(allNeuronAve_array_OppositeRF)
    set(gca,'XTick',[1 5 10]);
    set(gca, 'XTickLabels', {'-200', '0', '+250'});
    set(gca,'YTick', [4 8 12 16 20 24])
    set(gca, 'YTickLabels', {'76', '63', '50', '36', '23', '10'});
    title(s4, 'ArrayOnset - OppositeRF');
    ylabel('Frequence (Hz)');
    xlabel('Peri-array time (ms)')
    
    suptitle('Spike-field Coherence (all-neuron avg)');
    set(f1, 'Position', [452 136 1053 761]);
    cd(savedirOne)
    print(f1,'all-neuron avg', '-dpng');

    %    close;
    
    
    
    




