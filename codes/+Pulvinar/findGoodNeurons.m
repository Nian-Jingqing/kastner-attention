%% add dataset path
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
datasetPath = ['/snel/share/share/data/Kastner_Pulvinar/singleSession/continuous/data_raw/' ...
    'M20170608_PUL_all-g2-g3-g4-evokedSpiking-v8'];

%% Find file names and sort the files by modification time

filePattern = fullfile(datasetPath, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

cells = struct2cell(theFiles); % put the struct array into cell array
sortVals = cells(3,:)'; % take out the modification date from theFiles
ix = (1:length(sortVals))'; % initialize a vector for indexing the file sequence
time = datetime(sortVals); % convert modification time from cell to datetime for sorting
timeTable = table(time,ix); % put time and index in a table
sorted = sortrows(timeTable,'time'); % sort index by modification time 
theFiles = theFiles(sorted.ix); % use the index to sort theFiles

%% load the task info, spike times, and receptive fields for the neurons

d.units(length(theFiles)).rfLoc = zeros(1,2);
d.units(length(theFiles)).spikeTimes = 0;
% Initialize a 6x2 matrix to store the receptive fields (first column)
% and opposite fields (second column)

for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(datasetPath, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  if k == 1
      [d.UE]= Utils.loadTrialInfo(fullFileName, 'UE');
%       rt_TD = rt(isHoldTrial);
%       cueLoc_TD = cueLoc(isHoldTrial);
%       window_TD = (windowTargetDim(1) + windowTargetDim(2))*1000;
%       window_CO = (windowCueOnset(1) + windowCueOnset(2))*1000;
%       window_AO = (windowArrayOnset(1) + windowArrayOnset(2))*1000;
      % for the first round in the loop, load the basic info(response time, cued location,window)
  end
  [d.units(k).rfLoc, d.units(k).spikeTimes] = Utils.loadNeuron(fullFileName, 'spikeTs', 'inRFLoc', 'exRFLoc');
%   % load spiking data from files
%   spikeCounts_TD{k} = arrayfun(@(x) TrainToCounts(x.times,window_TD), spikeTimes_TD{k}, 'UniformOutput', false);
%   spikeCounts_CO{k} = arrayfun(@(x) TrainToCounts(x.times,window_CO), spikeTimes_CO{k}, 'UniformOutput', false);
%   spikeCounts_AO{k} = arrayfun(@(x) TrainToCounts(x.times,window_AO), spikeTimes_AO{k}, 'UniformOutput', false);
%   
%   spikeCountsAllTrialsMat_TD{k} = cell2mat((spikeCounts_TD{k})');
%   spikeCountsAllTrialsMat_CO{k} = cell2mat((spikeCounts_CO{k})');
%   spikeCountsAllTrialsMat_AO{k} = cell2mat((spikeCounts_AO{k})');
  
end

%% convert spikeTimes into 1s and 0s
sessStartTime = d.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1);
% session start time - I set it to be the time when the monkey starts
% fixation in the first trial
sessEndTime = d.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end);
% session end time - I set it to be the time when the monkey releases the
% lever in the last trial
total_samples = floor(1000 * (sessEndTime - sessStartTime))+1;
stream = [];
stream.spikes = sparse(total_samples, numel(d.units) );
rawSampleRate = 1000;
for iunit = 1:numel(d.units)
    spks = d.units(iunit).spikeTimes;
    spksshort = spks(spks > sessStartTime & spks < sessEndTime);
    spksshort = (spksshort - sessStartTime)*rawSampleRate;
    if abs(max(spksshort)-total_samples) < 1e-03
        spksshort(end) = spksshort(end) - 0.1;
    end
    flooredTrain = unique(floor(spksshort));
    stream.spikes(flooredTrain + 1, iunit) = 1;
end


%% extract trial information
startInds = round((d.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue - sessStartTime)*rawSampleRate);
startInds(1) = 1;
stopInds = round((d.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice - sessStartTime)*rawSampleRate);
trialstruct = struct;
for itrial = 1:numel(startInds)
    trialstruct(itrial).rt = d.UE.rt(itrial);
    trialstruct(itrial).cueLoc = d.UE.cueLoc(itrial);
    trialstruct(itrial).isHoldTrial = d.UE.isHoldTrial(itrial);
    trialstruct(itrial).startTime = d.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(itrial) - sessStartTime;
    trialstruct(itrial).endTime = d.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(itrial) - sessStartTime;
    trialstruct(itrial).startInd = startInds(itrial);
    trialstruct(itrial).endInd = stopInds(itrial);
end


%% put spike train into a Continuous class
dtMS = 1;
C = Continuous.Continuous(stream, dtMS);

%% turn into a trialized (R) struct
r = R.Rstruct( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

%% get first 25% of trials and last 25% of trials
r_copy = r.copy();
nTrials = numel(r_copy.r);
nNeurons = size(r.r(1).spikes,1);
percentSelected = 0.25;
trialNum = floor(nTrials*percentSelected);
r_copyFirst = R.Rstruct(r_copy.r(1:trialNum));
r_copyLast = R.Rstruct(r_copy.r(end-trialNum+1:end));

%% get avg firing rate for each neuron for each trial
avgFiringRateEachTrialFirst = zeros(nNeurons, trialNum);
avgFiringRateEachTrialLast = zeros(nNeurons, trialNum);

for itrial = 1:trialNum
    nTimesRawFirst = size(r_copyFirst.r(itrial).spikes,2);
    nTimesRawLast = size(r_copyLast.r(itrial).spikes,2);
    avgFiringRateEachTrialFirst(:,itrial) = sum(r_copyFirst.r(itrial).spikes, 2)/(nTimesRawFirst/1000);
    avgFiringRateEachTrialLast(:,itrial) = sum(r_copyLast.r(itrial).spikes, 2)/(nTimesRawLast/1000);
end

%% get avg firing rate for each neuron across the first 25% of the trials and the last 25% of the trials
avgFiringRateFirst = sum(avgFiringRateEachTrialFirst,2)*(1/trialNum);
avgFiringRateLast = sum(avgFiringRateEachTrialLast,2)*(1/trialNum);

%% plotting the avg firing rate (first 25% vs last 25%)

f1 = figure;
scatter(avgFiringRateFirst, avgFiringRateLast,10, 'b', 'filled');
hold on
x = linspace(0,max(avgFiringRateFirst));
y = x;
plot(x,y,'r', 'DisplayName', 'y = x');
hold on 
x1 = linspace(0,max(avgFiringRateFirst));
y1 = (4/5)*x1;
y2 = (6/5)*x1;
plot(x1,y1,'g', 'DisplayName', 'y = 0.8x');
hold on
plot(x1,y2,'g', 'DisplayName', 'y = 1.2x');
xlabel('Avg Firing Rate first 25% of trials');
ylabel('Avg Firing Rate last 25% of trials');
legend('show')
% cd(savedirOne);
set(f1, 'Position', [200 100 1200 1000]);
title('whole trials')
% print(f1,'cueOnset', '-dpng');

%% select good neurons

neuronIx = false(1,nNeurons);
for n = 1:nNeurons
    if (avgFiringRateLast(n) > avgFiringRateFirst(n)*(4/5)) && (avgFiringRateLast(n) < avgFiringRateFirst(n)*(6/5)) && avgFiringRateFirst(n) > 0.5
        neuronIx(n) = true;
    end
end


%%

%% add dataset path
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
datasetPath = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/datasets'];

savedirOne = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/findGoodNeurons/cueOnset'];
%% load dataset
datasetName = 'cueOnArrayOnTargetDim_HoldRel.mat';
fullFileName = fullfile(datasetPath, datasetName);
r_real = load(fullFileName); % get the original dataset (for all neurons)
r_real = R.Rstruct(r_real.R); % put the dataset into R struct class

%% get experiment info (nTrials, nTimes, nNeurons)
nTrials = length(r_real.r); % get trial number
nTimesRaw = size(r_real.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_real.r(1).spikeCounts, 1); % get neuron nubmer

%% get trials with the same alignment

r_realCopy = r_real.copy();
alignType = 2; % select arrayOnset type
alignIx = [r_realCopy.r.alignType];
r_realCopy = R.Rstruct(r_realCopy.r(alignIx == alignType));
nTrials = length(r_realCopy.r);
%% Select first 25% of the trials and last 25% of the trials

percentSelected = 0.25;
trialNum = floor(nTrials*percentSelected);
r_realCopyFirst = R.Rstruct(r_realCopy.r(1:trialNum));
r_realCopyLast = R.Rstruct(r_realCopy.r(end-trialNum+1:end));

%% Put the spike counts into matrix (nTrials x nTimes x nNeurons)

rawSpikingMatrixFirst = permute(cat(3, r_realCopyFirst.r.spikeCounts), [3 2 1]);
rawSpikingMatrixLast = permute(cat(3, r_realCopyLast.r.spikeCounts), [3 2 1]);

%% sum the spiking for all times and all trials for each neuron and get avg firing rate

avgFiringRateFirst = squeeze(sum(sum(rawSpikingMatrixFirst,1), 2)/(trialNum*(nTimesRaw/1000)));
avgFiringRateLast = squeeze(sum(sum(rawSpikingMatrixLast,1), 2)/(trialNum*(nTimesRaw/1000)));
avgFiringRate = [avgFiringRateFirst avgFiringRateLast];

%% plot the avg firing rate for the first 25% trials vs the avg firing rate for the last 25% trials for each neuron

f1 = figure;
scatter(avgFiringRate(:,1), avgFiringRate(:,2),10, 'b', 'filled');
hold on
x = linspace(0,max(avgFiringRate(:,1)));
y = x;
plot(x,y,'r', 'DisplayName', 'y = x');
hold on 
x1 = linspace(0,max(avgFiringRate(:,1)));
y1 = (4/5)*x1;
y2 = (6/5)*x1;
plot(x1,y1,'g', 'DisplayName', 'y = 0.8x');
hold on
plot(x1,y2,'g', 'DisplayName', 'y = 1.2x');
xlabel('Avg Firing Rate first 25% of trials');
ylabel('Avg Firing Rate last 25% of trials');
legend('show')
cd(savedirOne);
set(f1, 'Position', [200 100 1200 1000]);
title('alignedToCueOnset')
print(f1,'cueOnset', '-dpng');



%% select good neurons

neuronIx = false(1,nNeurons);
for n = 1:nNeurons
    if (avgFiringRate(n,2) > avgFiringRate(n,1)*(4/5)) && (avgFiringRate(n,2) < avgFiringRate(n,1)*(6/5)) && avgFiringRate(n,1) > 0.5
        neuronIx(n) = true;
    end
end

%% Use neuronIx to keep good neurons and save them to Rstruct





