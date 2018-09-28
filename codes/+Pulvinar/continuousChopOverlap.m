%% add dataset path
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')

%% find all the files to be processed
% Specify the folder where the files live.
dayName = 170621;
dayStr = num2str(dayName);
day_id = 7;
subFolderName = ['M20', dayStr];
rawDataBasePath = '/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/fourLocations/';
myFolder = fullfile(rawDataBasePath, subFolderName, 'MUA_GRATINGS');
%Specify where to save the reorganized data
saveDir = '/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/fourLocations/';
savefigDir = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
        'multi-unit/continuousOverlapChop/multiDay_JanToMar/externalInputsExamples_updated/170127/'];
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

theFiles = theFiles(nIndices);

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

%% load the task info, spike times, and receptive fields for the neurons

rfLoc = zeros(length(theFiles),2);
d.units(length(theFiles)).spikeTimes = 0;
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  
  if k == 1
      d.UE = loadInfoForContinuous(fullFileName);
  end
  
  [rfLoc(k,:), d.units(k).spikeTimes] = loadNeuronForContinuous(fullFileName);
  
end



%% convert spikeTimes into 1s and 0s
rawSampleRate = 1000;
sessStartTime = d.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1)*rawSampleRate;
% session start time - I set it to be the time when the monkey starts
% fixation in the first trial
sessEndTime = d.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
extraEndMs = 500;
%keep extra 500ms after each trial ends (level release)

% session end time - I set it to be the time when the monkey releases the
% lever in the last trial
total_samples = floor(sessEndTime - sessStartTime + extraEndMs)+1;
stream = [];
stream.spikes = sparse(total_samples, numel(d.units) );
for iunit = 1:numel(d.units)
    spks = d.units(iunit).spikeTimes;
    spksshort = spks(spks > sessStartTime & spks < sessEndTime);
    spksshort = spksshort - sessStartTime;
    if abs(max(spksshort)-total_samples) < 1e-03
        spksshort(end) = spksshort(end) - 0.1;
    end
    flooredTrain = unique(floor(spksshort));
    stream.spikes(flooredTrain + 1, iunit) = 1;
end

%% construct external input dimensions

%set cueonset external input timing
cueStart = round(d.UE.cueOnset * rawSampleRate - sessStartTime);
cueEnd = cueStart + 100;

%set arrayonset external input timing
arrayStart = round(d.UE.arrayOnset * rawSampleRate - sessStartTime);
arrayEnd = round(d.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice * rawSampleRate - sessStartTime);
arrayEnd_rel = arrayEnd(~d.UE.isHoldTrial);
%arrayEnd_hold = round(d.UE.targetDim * rawSampleRate + inputLength_hold - sessStartTime);
% arrayEnd = [arrayEnd_rel; arrayEnd_hold];
% arrayEnd = sort(arrayEnd);


%set targetDim external input timing
dimStart = round(d.UE.targetDimMatch(d.UE.isHoldTrial) * rawSampleRate - sessStartTime);
dimEnd = arrayEnd(d.UE.isHoldTrial);

% Figure out the array type based on isHoldTrial and cueLoc
arrayType = ones(length(arrayStart),1);
arrayPattern1Ix = (d.UE.isHoldTrial & (d.UE.cueLoc == 1 | d.UE.cueLoc == 3)) | (~d.UE.isHoldTrial & (d.UE.cueLoc == 2 | d.UE.cueLoc == 4));
arrayPattern2Ix = (d.UE.isHoldTrial & (d.UE.cueLoc == 2 | d.UE.cueLoc == 4)) | (~d.UE.isHoldTrial & (d.UE.cueLoc == 1 | d.UE.cueLoc == 3));
arrayType(arrayPattern2Ix) = 2;

%% initialize the externalInputs
cueLoc_list = unique(d.UE.cueLoc);
nCueLoc = length(cueLoc_list);
nExternalInputs = nCueLoc + 4;
stream.externalInputs = sparse(total_samples, nExternalInputs );  


%%
% set up cueOnset external inputs
for i = 1:nCueLoc
    cueStart_thisCue = cueStart(d.UE.cueLoc == cueLoc_list(i));
    cueEnd_thisCue = cueEnd(d.UE.cueLoc == cueLoc_list(i));
    for j = 1:length(cueStart_thisCue)
        stream.externalInputs(cueStart_thisCue(j) + 35:cueEnd_thisCue(j) + 35, i) = 1;
    end
end

% set up array external inputs for hold trials
arrayStart_hold = arrayStart(d.UE.isHoldTrial);
cueLoc_hold = d.UE.cueLoc(d.UE.isHoldTrial);
arrayChar_hold = d.UE.arrayShapesCorrect(d.UE.isHoldTrial);
for th = 1:length(arrayStart_hold)
    % starting at cueLoc, and follow the counter clockwise order for the
    % patterns
    stream.externalInputs(arrayStart_hold(th)+35:(dimStart(th)-1)+35, cueLoc_hold(th)+nCueLoc) = 1; % on cueLoc, must be hold shape
    stream.externalInputs(dimStart(th)+35:dimEnd(th)+35, cueLoc_hold(th)+nCueLoc) = 0.5; % Upon dimming, input intensity to 0.5
    
    if strcmp(arrayChar_hold(th), 'HHRR')
        % the 2nd character:
        stream.externalInputs(arrayStart_hold(th)+35:dimEnd(th)+35, cycleIndex(cueLoc_hold(th)+1)+nCueLoc) = 1;
        % the 3rd and 4th characters:
        stream.externalInputs(arrayStart_hold(th)+35:dimEnd(th)+35, [cycleIndex(cueLoc_hold(th)+2)+nCueLoc, cycleIndex(cueLoc_hold(th)+3)+nCueLoc]) = -1;
    end
    
    if strcmp(arrayChar_hold(th), 'HRHR')
        % the 3rd character:
        stream.externalInputs(arrayStart_hold(th)+35:dimEnd(th)+35, cycleIndex(cueLoc_hold(th)+2)+nCueLoc) = 1;
        % the 2nd and 4th characters:
        stream.externalInputs(arrayStart_hold(th)+35:dimEnd(th)+35, [cycleIndex(cueLoc_hold(th)+1)+nCueLoc, cycleIndex(cueLoc_hold(th)+3)+nCueLoc]) = -1;
    end
    
    if strcmp(arrayChar_hold(th), 'HRRH')
        % the 4th character:
        stream.externalInputs(arrayStart_hold(th)+35:dimEnd(th)+35, cycleIndex(cueLoc_hold(th)+3)+nCueLoc) = 1;
        % the 2nd and 3rd characters:
        stream.externalInputs(arrayStart_hold(th)+35:dimEnd(th)+35, [cycleIndex(cueLoc_hold(th)+1)+nCueLoc, cycleIndex(cueLoc_hold(th)+2)+nCueLoc]) = -1;
    end
end

% set up array external inputs for release trials
arrayStart_rel = arrayStart(~d.UE.isHoldTrial);
cueLoc_rel = d.UE.cueLoc(~d.UE.isHoldTrial);
arrayChar_rel = d.UE.arrayShapesCorrect(~d.UE.isHoldTrial);
for tr = 1:length(arrayStart_rel)
    % starting at cueLoc, and follow the counter clockwise order for the
    % patterns
    if strcmp(arrayChar_rel(tr), 'RRHH')
        % the 1st and 2nd character:
        stream.externalInputs(arrayStart_rel(tr)+35:arrayEnd_rel(tr)+35, [cueLoc_rel(tr)+nCueLoc, cycleIndex(cueLoc_rel(tr)+1)+nCueLoc]) = -1;
        % the 3rd and 4th characters:
        stream.externalInputs(arrayStart_rel(tr)+35:arrayEnd_rel(tr)+35, [cycleIndex(cueLoc_rel(tr)+2)+nCueLoc, cycleIndex(cueLoc_rel(tr)+3)+nCueLoc]) = 1;
    end
    
    if strcmp(arrayChar_rel(tr), 'RHRH')
        % the 1st and 3rd character:
        stream.externalInputs(arrayStart_rel(tr)+35:arrayEnd_rel(tr)+35, [cueLoc_rel(tr)+nCueLoc, cycleIndex(cueLoc_rel(tr)+2)+nCueLoc]) = -1;
        % the 2nd and 4th characters:
        stream.externalInputs(arrayStart_rel(tr)+35:arrayEnd_rel(tr)+35, [cycleIndex(cueLoc_rel(tr)+1)+nCueLoc, cycleIndex(cueLoc_rel(tr)+3)+nCueLoc]) = 1;
    end
    
    if strcmp(arrayChar_rel(tr), 'RHHR')
        % the 1st and 4th character:
        stream.externalInputs(arrayStart_rel(tr)+35:arrayEnd_rel(tr)+35, [cueLoc_rel(tr)+nCueLoc, cycleIndex(cueLoc_rel(tr)+3)+nCueLoc]) = -1;
        % the 2nd and 3rd characters:
        stream.externalInputs(arrayStart_rel(tr)+35:arrayEnd_rel(tr)+35, [cycleIndex(cueLoc_rel(tr)+1)+nCueLoc, cycleIndex(cueLoc_rel(tr)+2)+nCueLoc]) = 1;
    end
end

        
        

% % set up arrayOnset external inputs
% for i = 1:2
%     arrayStart_thisType = arrayStart(arrayType == i);
%     arrayEnd_thisType = arrayEnd(arrayType == i);
%     for j = 1:length(arrayStart_thisType)
%         stream.externalInputs(arrayStart_thisType(j):arrayEnd_thisType(j), (i + nCueLoc)) = 1;
%     end
% end
% 
% % set up array targetDim external inputs
% cueLoc_hold = d.UE.cueLoc(d.UE.isHoldTrial);
% for i = 1:nCueLoc
%     dimStart_thisCue = dimStart(cueLoc_hold == cueLoc_list(i));
%     dimEnd_thisCue = dimEnd(cueLoc_hold == cueLoc_list(i));
%     for j = 1:length(dimStart_thisCue)
%         stream.externalInputs(dimStart_thisCue(j):dimEnd_thisCue(j), (i + nCueLoc + 2)) = 1;
%     end
% end




%% extract trial information
startInds = round(d.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue * rawSampleRate - sessStartTime);
startInds(1) = 1;
stopInds = round(d.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice * rawSampleRate - sessStartTime)+400;
trialstruct = struct;
for itrial = 1:numel(startInds)
    trialstruct(itrial).rt = d.UE.rt(itrial);
    trialstruct(itrial).rfloc = rfLoc;
    trialstruct(itrial).cueLoc = d.UE.cueLoc(itrial);
    trialstruct(itrial).isHoldTrial = d.UE.isHoldTrial(itrial);
    trialstruct(itrial).isHoldBal = d.UE.isHoldBal(itrial);
    trialstruct(itrial).startTime = d.UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(itrial) * rawSampleRate - sessStartTime;
    trialstruct(itrial).endTime = d.UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(itrial) * rawSampleRate - sessStartTime;
    trialstruct(itrial).startInd = startInds(itrial);
    trialstruct(itrial).endInd = stopInds(itrial);
    
    % create condition type. This code is for 2 cueLoc. If using 4 cueLoc
    % data, then need to add stuff
    %if d.UE.cueLoc(itrial) == 1
    %   if arrayType == 1
    %       trialstruct(itrial).condition = 1;
    %   else
    %       trialstruct(itrial).condition = 2;
    %   end
    %elseif d.UE.cueLoc(itrial) == 3
    %   if arrayType == 1
    %       trialstruct(itrial).condition = 3;
    %   else
    %       trialstruct(itrial).condition = 4;
    %   end
    %end
    if d.UE.isHoldTrial(itrial)
        trialstruct(itrial).condition = d.UE.cueLoc(itrial);
    else
        trialstruct(itrial).condition = d.UE.cueLoc(itrial) + nCueLoc;
    end        
end


%% put spike train into a Continuous class
dtMS = 1;
C = Continuous.Continuous(stream, dtMS);
sigma_neural = 50;
C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );

%% turn into a trialized (R) struct
r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

%% Verify External inputs
% trial_selected = [1, 50, 100, 150, 200];
% dimLabel = round(d.UE.targetDimMatch * rawSampleRate - sessStartTime);
% 
% for i = 1: length(trial_selected)
%    
% % for i = 1
%     itrial = trial_selected(i);
%     ei = r.r(itrial).externalInputs;
%     x_cue = cueStart(itrial)-startInds(itrial);
%     x_cueEnd = cueEnd(itrial)-startInds(itrial);
%     x_array = arrayStart(itrial)-startInds(itrial);
%     x_dim = dimLabel(itrial)-startInds(itrial);
%     x_arrayEnd = arrayEnd(itrial)-1-startInds(itrial);
%     
%     f1 = figure;
%     
%     cue_1 = subplot(6,1,1);
%     plot(ei(1,:), 'r', 'DisplayName','Location 1 - cue');
%     if d.UE.isHoldTrial(itrial)
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_dim x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'targetDim', 'response'});
%         set(gca, 'FontSize', 7);
%        % title(['Trial ' int2str(itrial) ' - Hold'])
%     else
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'response'});
%         set(gca, 'FontSize', 7);
%        % title(['Trial ' int2str(itrial) ' - Release'])
%     end
%     ylim([-1 1]);
%     legend('show')
%     
%     cue_3 = subplot(6,1,2);
%     plot(ei(2,:), 'b', 'DisplayName','Location 3 - cue');
%     if d.UE.isHoldTrial(itrial)
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_dim x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'targetDim', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Hold'])
%     else
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'response'});
%         set(gca, 'FontSize', 7);
%        % title(['Trial ' int2str(itrial) ' - Release'])
%     end
%     ylim([-1 1]);
%     legend('show')
%     
%     array_1 = subplot(6,1,3);
%     plot(ei(3,:), 'r', 'DisplayName','Location 1 - array');
%     if d.UE.isHoldTrial(itrial)
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_dim x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'targetDim', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Hold'])
%     else
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Release'])
%     end
%     ylim([-1 1]);
%     legend('show')
%     
%     array_2 = subplot(6,1,4);
%     plot(ei(4,:), 'g', 'DisplayName','Location 2 - array');
%     if d.UE.isHoldTrial(itrial)
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_dim x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'targetDim', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Hold'])
%     else
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Release'])
%     end
%     ylim([-1 1]);
%     legend('show')
%     
%     array_3 = subplot(6,1,5);
%     plot(ei(5,:), 'b', 'DisplayName','Location 3 - array');
%     if d.UE.isHoldTrial(itrial)
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_dim x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'targetDim', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Hold'])
%     else
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Release'])
%     end
%     ylim([-1 1]);
%     legend('show')
%     
%     array_4 = subplot(6,1,6);
%     plot(ei(6,:), 'm', 'DisplayName','Location 4 - array');
%     if d.UE.isHoldTrial(itrial)
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_dim x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'targetDim', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Hold'])
%     else
%         set(gca,'XTick',[x_cue x_cueEnd x_array x_arrayEnd]);
%         set(gca,'XTickLabels',{'cue','+100ms','array', 'response'});
%         set(gca, 'FontSize', 7);
%         %title(['Trial ' int2str(itrial) ' - Release'])
%     end
%     ylim([-1 1]);
%     legend('show')
%     
%     if d.UE.isHoldTrial(itrial)
%         suptitle(['Trial ' int2str(itrial) ' - Hold'])
%     else
%         suptitle(['Trial ' int2str(itrial) ' - Release'])
%     end
%     
%     set(f1, 'Position', [135 86 1645 861]);
%     cd(savefigDir);
%     print(f1,['Trial ' int2str(itrial)], '-dpng');
%     close;
%     
% end
    
    

%% load pre-aligned data (for run manager to use the pre-aligned data to compute alignment Matrix)
preAlignedPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/preAligned/multi-day_CoAoTdHoldRel_JanToApr/withGoodNeurons_HoldRelSepForAO_fourLoc_v12'];
preAlignedDataFileName = [dayStr, '_cueOnArrayOnTargetDim_HoldRel.mat'];
preAlignedDir = fullfile(preAlignedPath, preAlignedDataFileName);

preAlign = load(preAlignedDir);

%%

combinedData.R = preAlign.R;
combinedData.r = r;
%% save
cd(saveDir);
saveName = [dayStr, '_cueOnArrayOnTargetDim_HoldRel.mat'];
save(saveName, 'combinedData');



% %%
% trial_time_ms = 500;
% trial_olap_ms = 100;
% data = r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );








