%% setup dates

%dataset(1).date = '02082019';
%dataset(2).date = '02132019';
%dataset(3).date = '02142019';
%dataset(4).date = '02152019';
%dataset(5).date = '02162019';
%dataset(6).date = '02262019';
%dataset(7).date = '02282019';
%dataset(8).date = '03012019';
%dataset(9).date = '03022019';
%dataset(10).date = '03032019';

%dataset(1).date = '02182019';
%dataset(2).date = '03062019';
%dataset(3).date = '03112019';
%dataset(4).date = '03142019';
%dataset(5).date = '03272019';
%dataset(6).date = '04062019';
%dataset(7).date = '04252019';
%dataset(8).date = '05022019';

%dataset(1).date = '02272019';
%dataset(2).date = '03042019';
%dataset(3).date = '03072019';
%dataset(4).date = '03092019';
%dataset(5).date = '03102019';
%dataset(6).date = '03122019';
%dataset(7).date = '03132019';
%dataset(8).date = '03152019';

dataset(1).date = '03162019';
dataset(2).date = '03182019';
dataset(3).date = '03292019';
dataset(4).date = '03312019';
dataset(5).date = '04012019';
dataset(6).date = '04032019';
dataset(7).date = '04052019';
dataset(8).date = '04242019';
dataset(9).date = '04262019';
dataset(10).date = '04292019';

%%
eventMatPath = '/snel/share/share/data/kastner/Manoj/LIP/eventMat/';
rawSampleRate = 1000;
savedir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/';

for day = 1:10
    tic;
    clear UE
    dayEventFileName = ['eventmat_' dataset(day).date '.mat'];
    fullEventFileName = fullfile(eventMatPath, dayEventFileName);
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
    saveFileName = ['UE_' dataset(day).date '.mat'];
    cd(savedir)
    save(saveFileName, 'UE');
    toc;
end
