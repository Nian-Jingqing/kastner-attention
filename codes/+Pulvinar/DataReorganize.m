%% add path

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

%% find all the files to be processed
% Specify the folder where the files live.
dayName = 170621;
dayStr = num2str(dayName);
day_id = 7;
subFolderName = ['M20', dayStr];
rawDataBasePath = '/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/fourLocations/';
myFolder = fullfile(rawDataBasePath, subFolderName, 'MUA_GRATINGS');
%Specify where to save the reorganized data
saveDir = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/preAligned/multi-day_CoAoTdHoldRel_JanToApr/withGoodNeurons_HoldRelSepForAO_fourLoc_v12'];

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
%%
fileInCell = struct2cell(theFiles);
fileNamesInCell = fileInCell(1, :);
[~, stringLength] = sort(cellfun(@length, fileNamesInCell), 'ascend');
theFiles = theFiles(stringLength);

% str=theFiles(1).name;
% tmp = strfind(str, 'M');
% cutStr = str(tmp(1):tmp(end));
% cutStr = str(tmp(1)+1:tmp(end)-1);
% tmp2=strfind(cutStr,'_');
% cutStr(tmp2(end)+1:end);
% str2num( cutStr(tmp2(end)+1:end) );


%%

% %% sort theFIles in the sequence based on modification time
% cells = struct2cell(theFiles); % put the struct array into cell array
% sortVals = cells(3,:)'; % take out the modification date from theFiles
% ix = (1:length(sortVals))'; % initialize a vector for indexing the file sequence
% time = datetime(sortVals); % convert modification time from cell to datetime for sorting
% timeTable = table(time,ix); % put time and index in a table
% sorted = sortrows(timeTable,'time'); % sort index by modification time 
% theFiles = theFiles(sorted.ix); % use the index to sort theFiles
% %% do this for 032917, keep the first half
% keep_indices = 1:2:127;
% theFiles = theFiles(keep_indices);
% 
% %% do this for 032917, keep the second half
% keep_indices = 2:2:128;
% theFiles = theFiles(keep_indices);

% %% do this for 033117, keep the smaller g
% keep_indices = 65:128;
% theFiles = theFiles(keep_indices);

% %% do this for 033117, keep the larger g
% keep_indices = 1:64;
% theFiles = theFiles(keep_indices);

%% Load the spike times and prefered and opposite locations for each neuron, and bin spike train to spike counts
% Also, load trial-specific info, such as response time, cue location, and
% window

rfLoc = zeros(length(theFiles),2);
% Initialize a 6x2 matrix to store the receptive fields (first column)
% and opposite fields (second column)

% initializing for targetDim data
spikeTimes_TD = cell(1,length(theFiles));
% Initialize a cell array to store the spike trains
spikeCounts_TD = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts
spikeCountsAllTrialsMat_TD = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts matrix for all trials
% for each neuron

% Initializing for cueOnset data
spikeTimes_CO = cell(1,length(theFiles));
% Initialize a cell array to store the spike trains
spikeCounts_CO = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts
spikeCountsAllTrialsMat_CO = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts matrix for all trials
% for each neuron

% Initializing for arrayOnset data
spikeTimes_AO = cell(1,length(theFiles));
% Initialize a cell array to store the spike trains
spikeCounts_AO = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts
spikeCountsAllTrialsMat_AO = cell(1,length(theFiles));
% Initialize a cell array to store the spike counts matrix for all trials
% for each neuron




for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  if k == 1
      [rt, cueLoc, isHoldTrial, isHoldBal, windowTargetDim, windowCueOnset, windowArrayOnset] = DataPreprocessing.loadInfo(fullFileName);
      rt_TD = rt(isHoldBal);
      cueLoc_TD = cueLoc(isHoldBal);
      window_TD = (windowTargetDim(1) + windowTargetDim(2))*1000;
      window_CO = (windowCueOnset(1) + windowCueOnset(2))*1000;
      window_AO = (windowArrayOnset(1) + windowArrayOnset(2))*1000;
      % for the first round in the loop, load the basic info(response time, cued location,window)
      
%       % for days 170329 and 170331
%       baseFileName_2 = theFiles(k+1).name;
%       fullFileName_2 = fullfile(myFolder, baseFileName_2);
%       [rt_2, cueLoc_2, isHoldTrial_2, windowTargetDim_2, windowCueOnset_2, windowArrayOnset_2] = DataPreprocessing.loadInfo(fullFileName_2);
%       rt_TD_2 = rt_2(isHoldTrial_2);
%       cueLoc_TD_2 = cueLoc_2(isHoldTrial_2);
%       rt = [rt; rt_2];
%       cueLoc = [cueLoc; cueLoc_2];
%       isHoldTrial = [isHoldTrial; isHoldTrial_2];
%       rt_TD = [rt_TD; rt_TD_2];
%       cueLoc_TD = [cueLoc_TD; cueLoc_TD_2];
      
     
  end
  
  [spikeTimes_TD{k}, spikeTimes_CO{k}, spikeTimes_AO{k}, rfLoc(k,:)]= DataPreprocessing.loadNeuron(fullFileName);
  % load spiking data from files
  spikeCounts_TD{k} = arrayfun(@(x) DataPreprocessing.TrainToCounts(x.times,window_TD), spikeTimes_TD{k}, 'UniformOutput', false);
  spikeCounts_CO{k} = arrayfun(@(x) DataPreprocessing.TrainToCounts(x.times,window_CO), spikeTimes_CO{k}, 'UniformOutput', false);
  spikeCounts_AO{k} = arrayfun(@(x) DataPreprocessing.TrainToCounts(x.times,window_AO), spikeTimes_AO{k}, 'UniformOutput', false);
  
  spikeCountsAllTrialsMat_TD{k} = cell2mat((spikeCounts_TD{k})');
  spikeCountsAllTrialsMat_CO{k} = cell2mat((spikeCounts_CO{k})');
  spikeCountsAllTrialsMat_AO{k} = cell2mat((spikeCounts_AO{k})');
  
end
% %% for 032917 and 033117, combine the two portions of trials together
% combinedSpikeCountsAllTrialsMat_TD = cell(1,0.5*length(theFiles));
% combinedSpikeCountsAllTrialsMat_CO = cell(1,0.5*length(theFiles));
% combinedSpikeCountsAllTrialsMat_AO = cell(1,0.5*length(theFiles));
% for i = 1:64
%     combinedSpikeCountsAllTrialsMat_TD{i} = cat(1, spikeCountsAllTrialsMat_TD{2*i - 1}, spikeCountsAllTrialsMat_TD{2*i});
%     combinedSpikeCountsAllTrialsMat_CO{i} = cat(1, spikeCountsAllTrialsMat_CO{2*i - 1}, spikeCountsAllTrialsMat_CO{2*i});
%     combinedSpikeCountsAllTrialsMat_AO{i} = cat(1, spikeCountsAllTrialsMat_AO{2*i - 1}, spikeCountsAllTrialsMat_AO{2*i});
% end
% TrialxTimexNeuron_TD = cat(3,combinedSpikeCountsAllTrialsMat_TD{:});
% TrialxTimexNeuron_CO = cat(3,combinedSpikeCountsAllTrialsMat_CO{:});
% TrialxTimexNeuron_AO = cat(3,combinedSpikeCountsAllTrialsMat_AO{:});

%% for other days
TrialxTimexNeuron_TD = cat(3,spikeCountsAllTrialsMat_TD{:});
TrialxTimexNeuron_CO = cat(3,spikeCountsAllTrialsMat_CO{:});
TrialxTimexNeuron_AO = cat(3,spikeCountsAllTrialsMat_AO{:});

%% Figure out the array type based on isHoldTrial and cueLoc
arrayType = ones(size(TrialxTimexNeuron_CO,1),1);
arrayPattern1Ix = (isHoldTrial & (cueLoc == 1 | cueLoc == 3)) | (~isHoldTrial & (cueLoc == 2 | cueLoc == 4));
arrayPattern2Ix = (isHoldTrial & (cueLoc == 2 | cueLoc == 4)) | (~isHoldTrial & (cueLoc == 1 | cueLoc == 3));
arrayType(arrayPattern2Ix) = 2;
arrayType_TD = arrayType(isHoldBal);



%% construct the R struct array
COCondStart = 1;
AOCondStart = length(unique(cueLoc))*1 + 1;
TDCondStart = length(unique(cueLoc))*3 + 1;
% length(unique(cueLoc)) tells you how many cue locations you have.
% Then, muliply number of cueloc by 2 (1 for CO and 1 for AO) tells you the
% stop position of the condition number for all trials for CO and AO

TD = struct;
for trial = 1:size(TrialxTimexNeuron_TD,1)
    TD(trial).spikeCounts = reshape(TrialxTimexNeuron_TD(trial,:,:),[size(TrialxTimexNeuron_TD,2),size(TrialxTimexNeuron_TD,3)])';
    TD(trial).rfloc = rfLoc;
    TD(trial).window = windowTargetDim;
    TD(trial).rt = rt_TD(trial);
    TD(trial).cueLoc = cueLoc_TD(trial);
    TD(trial).alignType = 3;
    TD(trial).arrayType = arrayType_TD(trial);
    TD(trial).isHoldTrial = 1;
    TD(trial).condition = TDCondStart + cueLoc_TD(trial)-1;

    %################## for 2-location data
    % if cueLoc_TD(trial) == 1
    %  TD(trial).condition = TDCondStart;
    %elseif cueLoc_TD(trial) == 3
    %       TD(trial).condition = TDCondStart + 1;
    %end
    %################## end
        
end

CO = struct;
for trial = 1:size(TrialxTimexNeuron_CO,1)
    CO(trial).spikeCounts = reshape(TrialxTimexNeuron_CO(trial,:,:),[size(TrialxTimexNeuron_CO,2),size(TrialxTimexNeuron_CO,3)])';
    CO(trial).rfloc = rfLoc;
    CO(trial).window = windowCueOnset;
    CO(trial).rt = rt(trial);
    CO(trial).cueLoc = cueLoc(trial);
    CO(trial).alignType = 1;
    CO(trial).arrayType = arrayType(trial);
    CO(trial).isHoldTrial = isHoldTrial(trial);
    CO(trial).condition = COCondStart + cueLoc(trial)-1;
    %################ for 2-location data
    %if cueLoc(trial) == 1
    %   CO(trial).condition = COCondStart;
    %elseif cueLoc(trial) == 3
    %       CO(trial).condition = COCondStart + 1;
    %end
end

AO = struct;
for trial = 1:size(TrialxTimexNeuron_AO,1)
    AO(trial).spikeCounts = reshape(TrialxTimexNeuron_AO(trial,:,:),[size(TrialxTimexNeuron_AO,2),size(TrialxTimexNeuron_AO,3)])';
    AO(trial).rfloc = rfLoc;
    AO(trial).window = windowArrayOnset;
    AO(trial).rt = rt(trial);
    AO(trial).cueLoc = cueLoc(trial);
    AO(trial).alignType = 2;
    AO(trial).arrayType = arrayType(trial);
    AO(trial).isHoldTrial = isHoldTrial(trial);
    if isHoldTrial(trial)
        AO(trial).condition = AOCondStart + cueLoc(trial)-1;
    else
        AO(trial).condition = AOCondStart + length(unique(cueLoc))+ cueLoc(trial)-1;
    end
    %################# for 2-location data ###########
    %if isHoldTrial(trial)    
    %   if cueLoc(trial) == 1
    %       AO(trial).condition = AOCondStart;
    %   elseif cueLoc(trial) == 3
    %       AO(trial).condition = AOCondStart + 1;
    %   end
    %elseif ~isHoldTrial(trial)
    %   if cueLoc(trial) == 1
    %       AO(trial).condition = AOCondStart + 2;
    %   elseif cueLoc(trial) == 3
    %       AO(trial).condition = AOCondStart + 3;
    %   end
    % end
end

R = [CO,AO,TD];

% %% select good neurons - runs until 03/02/2018
% 
% % nIndices = [1 3 4 6 7 8 11 12 13 14 15 16 17 19 20 21 22 23 24 25 26 28 29 30 31 32];
% % 170127
% % nIndices = [1 2 4 5 6 7 8 9 10 11 12 14 15 16 17 19 20 21 23 24 26 29 30 31 32];
% % % 170130
% % nIndices = [1 5 8 10 12 13 21 26 27 28 30 32];
% % % 170201
% % nIndices = [3 4 5 7 11 12 15 16 17 18 23 28 30];
% % % 170211
% % nIndices = [9 15 18 19 20 21 26 28 30 31 32 33 34 35 37 38 39 41 42 46 49 55 58 60 62 63 64];
% % % 170308
% % nIndices = [1 2 3 5 8 9 11 13 15 18 20 21 22 23 32 35 38 42 43 45 46 47 48 50 54 55 56 57 58 59 63];
% % % 170311
% for i = 1:length(R)
%     R(i).spikeCounts = R(i).spikeCounts(nIndices,:);
%     R(i).rfloc = R(i).rfloc(nIndices,:);
% end


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
% For all four targets - v12 (for runs after 07/17/2018)
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


%% select good neurons v10 - 03/12/2018 - 05/21/18

% nIndices = [1 3 4 6 7 8 11 13 14 15 16 17 19 20 21 22 23 24 25 26 29 30 31 32];
% % 170127
% nIndices = [1 2 5 6 7 8 9 12 14 15 16 17 19 20 21 23 24 26 29 30 31];
% % 170130
% nIndices = [1 5 8 10 12 13 27 30 32];
% % 170201
% nIndices = [3 4 5 7 11 12 15 16 17 18 23 27 28 30];
% % 170211
% nIndices = [9 15 18 19 20 28 30 31 32 33 34 37 38 39 42 46 49 55 58 62 63 64];
% % 170308
% nIndices = [1 2 3 8 9 11 13 18 20 22 32 35 38 42 43 45 48 50 54 55 56 59];
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
%%
for i = 1:length(R)
    R(i).spikeCounts = R(i).spikeCounts(nIndices,:);
    R(i).rfloc = R(i).rfloc(nIndices,:);
end

% %% chop the trials into desired length: in this case, chop to 400ms before to 400ms after the event
% 
% nTimes = window_CO;
% chopStart = 1*(nTimes/4) + 1; 
% chopEnd = nTimes - 1*(nTimes/4);
% for i = 1:length(R)
%     R(i).spikeCounts = R(i).spikeCounts(:, chopStart:chopEnd);
%     R(i).window = [0.4, 0.4];
% end


%% save
cd(saveDir);
saveName = [dayStr, '_cueOnArrayOnTargetDim_HoldRel.mat'];
save(saveName, 'R');

    
    
% %% rand suffle neurons and make held-out neurons for GLM
% R_all = R;
% rng('default');
% for i = 1:5
%     R = R_all;
%     x = randperm(106);
%     train = x(1:75);
%     test = x(76:end);
%     for unit = 1:length(R_all)
%         R(unit).spikeCounts = R_all(unit).spikeCounts(train,:);
%         R(unit).rfloc = R_all(unit).rfloc(train,:);
%     end
%     cd(saveDir);
%     formatSpec = 'cueOnArrayOnTargetDim_HoldRel_00%d.mat';
%     datasetName = sprintf(formatSpec,i);
%     save(datasetName, 'R');
% end
%     

% %% CP - sort by file name - 03/14/2018
% 
% str=theFiles(1).name;
% tmp = strfind(str, 'M');
% cutStr = str(tmp(1):tmp(end));
% cutStr = str(tmp(1)+1:tmp(end)-1);
% tmp2=strfind(cutStr,'_');
% cutStr(tmp2(end)+1:end);
% str2num( cutStr(tmp2(end)+1:end) );




