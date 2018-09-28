%% add dataset path
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')

%% set directory

myFolder = '/snel/share/share/data/kastner/pulvinar/lfp/M20170311/';
%Specify where to save the reorganized data
saveDir = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/continuousOverlapChop/170311/';


%% load data

%{
%fileName = 'M20170127-ch1-ch32-g3-g4-g5-continuousLfpsSpikes-v12.mat';
%fullDir = fullfile(myFolder, fileName);
%data = load(fullDir);
%}

%%
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFile = dir(filePattern);
fileName = theFile(2).name;
fullDir = fullfile(myFolder, fileName);
data = load(fullDir);
% %% load timestamps from previous data file
% 
% loadpath = ['/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/' ...
%     'M20170127/Gratings/M20170127_PUL_1M-g3-g4-g5-evokedSpiking-v10.mat'];
% prev_data = load(loadpath);
% UE_longerSession = prev_data.UE;
% clear prev_data

% %% calculate offset and substract the offset from the timestamps retrieved from previous data files
% 
% offset = UE_longerSession.cueOnset - data.UE.cueOnset;
% 
% 
% firstEnterFixationTimesPreCue = UE_longerSession.fixationAndLeverTimes.firstEnterFixationTimesPreCue - offset;
% firstLeverReleaseTimesAroundJuice = UE_longerSession.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice - offset;

%% find sessStartTime and sessEndTimes
rawSampleRate = 1000;
sessStartTime = data.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1)*rawSampleRate;
% session start time - I set it to be the time when the monkey starts
% fixation in the first trial
sessEndTime = data.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
extraEndMs = 500;

sessChopStart = floor(sessStartTime) + 1;
total_samples = floor(sessEndTime - sessStartTime + extraEndMs)+1;
stream = [];
stream.lfps = data.channelDataNorm(:, sessChopStart:sessChopStart + total_samples - 1);

% %% get rid of the common average channel
% stream.lfps = stream.lfps(1:32, :);

% %% bandpass LFP data in the continuous form
% filtHighCutoff = 15;
% filtLowCutoff = 2;
% Fs = 1000;
% stream.lfps = bandpassFilter_continuous(stream.lfps, filtHighCutoff, filtLowCutoff, Fs);
stream.lfps = stream.lfps'; % turn lfp into the format for continuous class

% ############ using offset #################################

% %% calculate indices for sessStartTimeand sessEndTime
% rawSampleRate = 1000;
% 
% sessStartTime_longerSession =  UE_longerSession.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1)*rawSampleRate;
% sessEndTime_longerSession = UE_longerSession.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
% 
% sessStartTime = firstEnterFixationTimesPreCue(1)*rawSampleRate;
% sessEndTime = firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
% extraEndMs = 500;
% 
% %% chop the LFP data to make the LFP data start from sessStartTime
% sessChopStart = floor(sessStartTime) + 1;
% total_samples = floor(sessEndTime_longerSession - sessStartTime_longerSession + extraEndMs)+1;
% % to make sure the size of total_samples for lfp matches that for spikes
% % that was processed
% stream = [];
% stream.lfps = data.channelDataNorm(:, sessChopStart:sessChopStart + total_samples - 1);
% 
%##############################################################



%% extract trial information
startInds = round(data.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue * rawSampleRate - sessStartTime);
startInds(1) = 1;
stopInds = round(data.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice * rawSampleRate - sessStartTime)+400;
trialstruct = struct;

holdtrialCount = 0;
for itrial = 1:numel(startInds)
    trialstruct(itrial).rt = data.UE.rt(itrial);
%     trialstruct(itrial).rfloc = rfLoc;
    trialstruct(itrial).cueLoc = data.UE.cueLoc(itrial);
    trialstruct(itrial).isHoldTrial = data.UE.isHoldTrial(itrial);
    trialstruct(itrial).isTrialOutlier = data.isTrialOutlier(itrial);
    trialstruct(itrial).isNoisyChannel = data.isNoisyChannel;
    trialstruct(itrial).startTime = data.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(itrial) * rawSampleRate - sessStartTime;
    trialstruct(itrial).endTime = data.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(itrial) * rawSampleRate - sessStartTime;
    trialstruct(itrial).cueOnsetTime = data.UE.cueOnset(itrial) * rawSampleRate - sessStartTime;
    trialstruct(itrial).arrayOnsetTime = data.UE.arrayOnset(itrial) * rawSampleRate - sessStartTime;
    if data.UE.isHoldTrial(itrial)
        holdtrialCount = holdtrialCount + 1;
        trialstruct(itrial).targetDimTime = data.UE.targetDim(holdtrialCount) * rawSampleRate - sessStartTime;
    else
        trialstruct(itrial).targetDimTime = nan;
    end
    trialstruct(itrial).startInd = startInds(itrial);
    trialstruct(itrial).endInd = stopInds(itrial);
    
%     % create condition type. This code is for 2 cueLoc. If using 4 cueLoc
%     % data, then need to add stuff
%     if d.UE.cueLoc(itrial) == 1
%         if arrayType == 1
%             trialstruct(itrial).condition = 1;
%         else
%             trialstruct(itrial).condition = 2;
%         end
%     elseif d.UE.cueLoc(itrial) == 3
%         if arrayType == 1
%             trialstruct(itrial).condition = 3;
%         else
%             trialstruct(itrial).condition = 4;
%         end
%     end
            
end


%% put spike train into a Continuous class
dtMS = 1;
C = Continuous.Continuous(stream, dtMS);
% sigma_neural = 50;
% C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );

%% turn into a trialized (R) struct
r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

%% load the first half of data
firstHalfDataDir = fullfile(saveDir, '170311_cueOnArrayOnTargetDim_HoldRel_lfp_firstHalf.mat');
firstHalfData = load(firstHalfDataDir);
r_firstHalf = firstHalfData.r;

%% combine the first half and second half
for i = 1: length(r.r)
    r.r(i).lfps = [r_firstHalf.r(i).lfps(1:end-1,:); r.r(i).lfps];
    r.r(i).isNoisyChannel = [r_firstHalf.r(i).isNoisyChannel(1:end-1,:); r.r(i).isNoisyChannel];
end

%% save
cd(saveDir);
save('170311_cueOnArrayOnTargetDim_HoldRel_lfp.mat', 'r');

    
    

