%% find all the files to be processed
% Specify the folder where the files live.
myFolder = '/Users/feng/SNEL/Computational/CollaboratingProject/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7';
%Specify where to save the reorganized data
saveDir = '/Users/feng/SNEL/Computational/CollaboratingProject/ReorganizedData/TargetDim';

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);


%% Load the spike times and prefered and opposite locations for each neuron, and bin spike train to spike counts
% Also, load trial-specific info, such as response time, cue location, and
% window
spikeTimes = cell(1,length(theFiles));
% Initialize a cell array to store the spike trains
spikeCounts = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts
spikeCountsAllTrialsMat = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts matrix for all trials
% for each neuron
rfLoc = zeros(length(theFiles),2);
% Initialize a 6x2 matrix to store the receptive fields (first column)
% and opposite fields (second column)


for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  if k == 1
      [rt, cueLoc, windowTargetDim]= loadInfo(fullFileName);
      window = (windowTargetDim(1) + windowTargetDim(2))*1000;
      % for the first round in the loop, load the basic info(response time, cued location,window)
  end
  [spikeTimes{k}, rfLoc(k,:)]= loadNeuron(fullFileName);
  spikeCounts{k} = arrayfun(@(x) TrainToCounts(x.times,window), spikeTimes{k}, 'UniformOutput', false);
  spikeCountsAllTrialsMat{k} = cell2mat((spikeCounts{k})');
end
TrialxTimexNeuron = cat(3,spikeCountsAllTrialsMat{:});

%% construct the R struct array

R = struct;
for trial = 1:size(TrialxTimexNeuron,1)
    R(trial).spikeCounts = reshape(TrialxTimexNeuron(trial,:,:),[size(TrialxTimexNeuron,2),size(TrialxTimexNeuron,3)])';
    R(trial).rfloc = rfLoc;
    R(trial).window = windowTargetDim;
    R(trial).rt = rt(trial);
    R(trial).cueLoc = cueLoc(trial);
end
cd(saveDir);
save('TargetDim_holdTrials.mat', 'R');
    
    
    