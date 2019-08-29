%% add dataset path
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')

%% set up day information
datasets(1).shortName = '0208';
datasets(1).longName = '02082019';
datasets(2).shortName = '0218';
datasets(2).longName = '02182019';
datasets(3).shortName = '0226';
datasets(3).longName = '02262019';
datasets(4).shortName = '0227';
datasets(4).longName = '02272019';
datasets(5).shortName = '0308';
datasets(5).longName = '03082019';
datasets(6).shortName = '0310';
datasets(6).longName = '03102019';
datasets(7).shortName = '0311';
datasets(7).longName = '03112019';

%%
% assign rawDataRoot
rawDataRoot = '/snel/share/share/data/kastner/Manoj/PUL/spikeSorted';

%for day = 1:numel(datasets)
for day = 2

    %% construct the data path for loading
    folderName = ['Remy_' datasets(day).shortName '_PUL'];
    fileName = ['Remy_' datasets(day).longName '_PUL_MUA.mat'];
    fullFileName = fullfile(rawDataRoot, folderName, fileName);

    %% load the data
    S = load(fullFileName);

    %% combine units on the same channel
    % get rid of other channels that are not spiking channels
    allFieldNames = fieldnames(S);
    unit_names = {};
    j = 1;
    for i = 1:numel(allFieldNames)
        if strcmp(allFieldNames{i}(1:3), 'SPK')
            unit_names{j} = allFieldNames{i};
            j = j+1;
        end
    end

    concat_allUnits = [];
    allChannels_spks = {};
    sortOrNot = 0;
    j = 1;
    for i = 1:numel(unit_names)
        if i == 1
            tmp_0 = strfind(unit_names{i}, 'C');
            curr_Str = unit_names{i}(tmp_0 + 1 : tmp_0 + 3);
            curr_unit = str2num(curr_Str);
            concat_allUnits = S.(unit_names{i});
        elseif i == numel(unit_names) & sortOrNot
            allChannels_spks{j} = sort(concat_allUnits);
            break
        elseif i == numel(unit_names) & ~sortOrNot
            allChannels_spks{j} = concat_allUnits;
            break
        end
        tmp = strfind(unit_names{i+1}, 'C');
        next_Str = unit_names{i+1}(tmp + 1 : tmp + 3);
        next_unit = str2num(next_Str);    
        if curr_unit == next_unit
            concat_allUnits = [concat_allUnits; S.(unit_names{i+1})];
            sortOrNot = 1;
        elseif sortOrNot
            allChannels_spks{j} = sort(concat_allUnits);
            sortOrNot = 0;
            j = j+1;
            concat_allUnits = S.(unit_names{i+1});
        else
            allChannels_spks{j} = concat_allUnits;
            sortOrNot = 0;
            j = j+1;
            concat_allUnits = S.(unit_names{i+1});
        end
        curr_unit = next_unit;
    end

    % get rid of channels 1 - 4
    allChannels_spks = allChannels_spks(3:end);

    %% convert spikeTimes into 1s and 0s
    rawSampleRate = 1000;
    % in this dataset, all eventtimes and spike times are aligned to Start, which is 0. 
    sessStartTime = S.Start*rawSampleRate;
    sessEndTime = S.Stop*rawSampleRate;

    %extraEndMs = 500;
    %keep extra 500ms after each trial ends (level release)

    %total_samples = floor(sessEndTime - sessStartTime + extraEndMs)+1;
    total_samples = floor(sessEndTime - sessStartTime)+1;
    stream = [];
    stream.spikes = sparse(total_samples, numel(allChannels_spks) );
    for iunit = 1:numel(allChannels_spks)
        spks = allChannels_spks{iunit}*rawSampleRate;
        spksshort = spks(spks > sessStartTime & spks < sessEndTime);
        spksshort = spksshort - sessStartTime;
        if abs(max(spksshort)-total_samples) < 1e-03
            spksshort(end) = spksshort(end) - 0.1;
        end
        flooredTrain = unique(floor(spksshort));
        stream.spikes(flooredTrain + 1, iunit) = 1;
    end

    %% load event times
    % dayEventFileName = 'eventmat.mat';
    %dayEventFileName = 'Remy_eventmat_0310.mat';
    dayEventFileName = ['Remy_eventmat_' datasets(day).shortName '.mat'];
    fullEventFileName = fullfile(rawDataRoot, folderName, dayEventFileName);
    %fullEventFileName = fullfile(rawDataPath, dayEventFileName);
    load(fullEventFileName);
    aet = eventmatT; % aet = all event times
    aet.times = aet.times * rawSampleRate;
    %%
    % get all event codes for targetOnset and put in a vector, but exclude saccade error trials
    targetCodes = [20, 21,22, 28, 29, 30, 54, 55, 56, 62, 63, 64];
    saccadeErrorCodes = [70, 72];
    %totalNumTrials = sum(sum(aet.code == targetCodes));
    targetOnIndices_all = find(ismember(aet.code, targetCodes));
    targetOnIndices = targetOnIndices_all(~ismember(aet.code(targetOnIndices_all+1), saccadeErrorCodes));
    totalNumTrials = numel(targetOnIndices);

    %% establish a UE struct to store trial information
    UE.isValidTarget = false(1,totalNumTrials);
    UE.isSameObjTarget = false(1,totalNumTrials);
    UE.isDiffObjTarget = false(1,totalNumTrials);

    UE.cueType = zeros(1,totalNumTrials);
    UE.barType = zeros(1,totalNumTrials);
    UE.fixType = zeros(1,totalNumTrials);

    UE.isErrorTrial = NaN(1,totalNumTrials);
    UE.isEarlyError = zeros(1,totalNumTrials); % for trials that the release within 100ms of target offset
                                               % this is different from the codes 32, 66 in the event codes. Manoj includes trials that release even before target onset in these codes 

    UE.fixOn = zeros(1,totalNumTrials);
    UE.barOn = zeros(1,totalNumTrials);
    UE.cueOn = zeros(1,totalNumTrials);
    UE.targetOn = zeros(1,totalNumTrials);
    UE.trialEnd = zeros(1,totalNumTrials);

    cueCodes = [16:19, 24:27, 50:53, 58:61];
    barCodes = [15, 23, 49, 57];
    fixCodes = [14, 48];
    outcomeCodes = [31:34, 65:68];

    for i = 1 : numel(targetOnIndices)
        UE.targetOn(i) = round(aet.times(targetOnIndices(i)));
        UE.isValidTarget(i) = ismember(aet.code(targetOnIndices(i)), [20, 28, 54, 62]);
        UE.isSameObjTarget(i) = ismember(aet.code(targetOnIndices(i)), [21, 29, 55, 63]);
        UE.isDiffObjTarget(i) = ismember(aet.code(targetOnIndices(i)), [22, 30, 56, 64]);
        if ismember(aet.code(targetOnIndices(i) - 1), cueCodes)
            UE.cueOn(i) = round(aet.times(targetOnIndices(i) - 1));
            if ismember(aet.code(targetOnIndices(i) - 1), [16,24]) % group hori and vert together, Exo TL
                UE.cueType(i) = 1;
            end
            if ismember(aet.code(targetOnIndices(i) - 1), [17,25]) % Exo BL
                UE.cueType(i) = 2;
            end
            if ismember(aet.code(targetOnIndices(i) - 1), [18,26]) % Exo TR
                UE.cueType(i) = 3;
            end
            if ismember(aet.code(targetOnIndices(i) - 1), [19,27]) % Exo BR
                UE.cueType(i) = 4;
            end
            if ismember(aet.code(targetOnIndices(i) - 1), [50,58]) % Endo TL
                UE.cueType(i) = 5;
            end
            if ismember(aet.code(targetOnIndices(i) - 1), [51,59]) % Endo BL
                UE.cueType(i) = 6;
            end
            if ismember(aet.code(targetOnIndices(i) - 1), [52,60]) % Endo TR
                UE.cueType(i) = 7;
            end
            if ismember(aet.code(targetOnIndices(i) - 1), [53,61]) % Endo BR
                UE.cueType(i) = 8;
            end        

        else
            disp('Not a cue code')
        end

        if ismember(aet.code(targetOnIndices(i) - 2), barCodes)
            UE.barOn(i) = round(aet.times(targetOnIndices(i) - 2));
            if ismember(aet.code(targetOnIndices(i) - 2), [15, 49]) % vertical bars
                UE.barType(i) = 1;
            elseif ismember(aet.code(targetOnIndices(i) - 2), [23, 57]) % horizontal bars
                UE.barType(i) = 2;
            end
        else
            disp('Not a bar code')
        end
        
        if ismember(aet.code(targetOnIndices(i) - 3), fixCodes)
            UE.fixOn(i) = round(aet.times(targetOnIndices(i) - 3));
            if aet.code(targetOnIndices(i) - 3) == 14 % exo
                UE.fixType(i) = 1;
            elseif aet.code(targetOnIndices(i) - 3) == 48 % endo
            UE.fixType(i) = 2;
            end        
        else
            disp('Not a fix code')
        end

        if ismember(aet.code(targetOnIndices(i) + 1), outcomeCodes)
            UE.trialEnd(i) = round(aet.times(targetOnIndices(i) + 1));
            if ismember(aet.code(targetOnIndices(i) + 1), [31, 33, 65, 67]) % no bar response combined no response at all and no response during the window
                UE.isErrorTrial(i) = 1;
            elseif ismember(aet.code(targetOnIndices(i) + 1), [32, 66]) % bar response within 100ms of the target Offset
                UE.isErrorTrial(i) = 1;
                UE.isEarlyError(i) = 1;
            elseif ismember(aet.code(targetOnIndices(i) + 1), [34, 68])
                UE.isErrorTrial(i) = 0;
            end        
        else
            disp('Not a outcome code')
        end
        
    end

    %% construct external inputs
    UE.cueEnd = UE.cueOn + 150;
    UE.targetEnd = UE.targetOn + 150;
    UE.targetEnd(UE.targetEnd >= UE.trialEnd) = UE.trialEnd(UE.targetEnd >= UE.trialEnd);
    % make sure the targetEnd is the same as the trialEnd for the early stopping trials

    nCueTypes = numel(unique(UE.cueType));
    nFixTypes = numel(unique(UE.fixType));
    nBarTypes = numel(unique(UE.barType));
    nTargetTypes = 4; 

    fixDims = sparse(total_samples, 1 + nFixTypes);
    barDims = sparse(total_samples, nBarTypes);
    cueDims = sparse(total_samples, nCueTypes/2); % because these cueDims are only for Exo
    targetDims = sparse(total_samples, nTargetTypes);
    % construct fix square dimensions
    % Total five dimensions. First dimension - plain color. The rest four dimensions - cue colors

    for i = 1:numel(targetOnIndices)
        
        % construct fix square dimensions
        % Total five dimensions. First dimension - plain color. The rest four dimensions - cue colors
        if UE.fixType(i) == 1
            fixDims(UE.fixOn(i) + 35 : UE.trialEnd(i) + 35, 1) = 1; % first dimension is plain color
        elseif UE.fixType(i) == 2 % for color cue, turn down first dim during the cue
            fixDims(UE.fixOn(i) + 35 : UE.cueOn(i) + 35, 1) = 1;
            fixDims(UE.cueEnd(i) + 35 : UE.trialEnd(i) + 35, 1) = 1;
            fixDims(UE.cueOn(i) + 35 : UE.cueEnd(i) + 35, UE.cueType(i) - (nCueTypes/2 - 1)) = 1;
        end

        % construct bar dimensions

        barDims(UE.barOn(i) + 35 : UE.trialEnd(i) + 35, UE.barType(i)) = 1;

        % construct cue dimensions (only for exo cues, because endo are changes in fix))
        if UE.fixType(i) == 1
            cueDims(UE.cueOn(i) + 35 : UE.cueEnd(i) + 35, UE.cueType(i)) = 1;
        end

        % construct target dimensions
        targetPos = DataPreprocessing.Manoj_getTargetPos(UE.cueType(i));
        if UE.isValidTarget(i)
            targetDims(UE.targetOn(i) + 35 : UE.targetEnd(i) + 35, targetPos) = 1;
        elseif (UE.isSameObjTarget(i) && UE.barType(i) == 1) || (UE.isDiffObjTarget(i) && UE.barType(i) == 2)
            if ismember(targetPos, [1, 3])
                shiftedTargetPos = targetPos + 1;
            else
                shiftedTargetPos = targetPos - 1;
            end
            targetDims(UE.targetOn(i) + 35 : UE.targetEnd(i) + 35, shiftedTargetPos) = 1;
        elseif (UE.isSameObjTarget(i) && UE.barType(i) == 2) || (UE.isDiffObjTarget(i) && UE.barType(i) == 1)
            if ismember(targetPos, [1, 2])
                shiftedTargetPos = targetPos + 2;
            else
                shiftedTargetPos = targetPos - 2;
            end        
            targetDims(UE.targetOn(i) + 35 : UE.targetEnd(i) + 35, shiftedTargetPos) = 1;
        end

    end
    %%
    stream.externalInputs = [fixDims, barDims, cueDims, targetDims];


    %% extract trial information
    startInds = UE.fixOn;
    stopInds = UE.trialEnd + 400;
    trialstruct = struct;
    for itrial = 1:numel(startInds)
        trialstruct(itrial).isErrorTrial = UE.isErrorTrial(itrial);
        trialstruct(itrial).isEarlyError = UE.isEarlyError(itrial);
        trialstruct(itrial).cueType = UE.cueType(itrial);
        trialstruct(itrial).barType = UE.cueType(itrial);
        trialstruct(itrial).fixType = UE.fixType(itrial);
        trialstruct(itrial).isValidTarget = UE.isValidTarget(itrial);
        trialstruct(itrial).isSameObjTarget = UE.isSameObjTarget(itrial);
        trialstruct(itrial).isDiffObjTarget = UE.isDiffObjTarget(itrial);    
        trialstruct(itrial).startTime = startInds(itrial);
        trialstruct(itrial).endTime = stopInds(itrial);
        trialstruct(itrial).startInd = startInds(itrial);
        trialstruct(itrial).endInd = stopInds(itrial);
        trialstruct(itrial).condition = UE.cueType(itrial); % need to change!!        
    end

    %% put spike train into a Continuous class
    dtMS = 1;
    C = Continuous.Continuous(stream, dtMS);
    sigma_neural = 10;
    C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );

    %% turn into a trialized (R) struct
    r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

    %% combine and save
    combinedData.r = r;
    combinedData.UE = UE;
    %% save
    saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/';
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    cd(saveDir);
    saveName = [datasets(day).shortName '19_v1.mat'];
    %saveName = '020819_v1.mat';
    save(saveName, 'combinedData');

end
