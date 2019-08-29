%% list the dates for all needed datasets
datasets(1).shortName = '0218';
datasets(1).midName = '021819';
datasets(1).longName = '02182019';
datasets(2).shortName = '0306';
datasets(2).midName = '030619';
datasets(2).longName = '03062019';
datasets(3).shortName = '0311';
datasets(3).midName = '031119';
datasets(3).longName = '03112019';
datasets(4).shortName = '0314';
datasets(4).midName = '031419';
datasets(4).longName = '03142019';
datasets(5).shortName = '0406';
datasets(5).midName = '040619';
datasets(5).longName = '04062019';
datasets(6).shortName = '0425';
datasets(6).midName = '042519';
datasets(6).longName = '04252019';
datasets(7).shortName = '0502';
datasets(7).midName = '050219';
datasets(7).longName = '05022019';

%% get drifting periods and bad channels due to drifting
drifty_info


%% bad channels due to high correlation
hc_channels{1} = [24, 26];
hc_channels{2} = [24, 25, 27];
hc_channels{3} = [9];
hc_channels{4} = [];
hc_channels{5} = [18];
hc_channels{6} = [];
hc_channels{7} = [32];
allChannels = 1:32;


%% % load all the chopped and recombined data (post-LFADS)
loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/withoutDataMasking/selected/';
% iterate over days, load each day and add it to a 'olapChopped' cell array
clear tc_origin
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_v1.mat', loadpath, datasets( nday ).midName );
    tmp = load( fname );
    tc_origin{ nday } = tmp.combinedData;
    %UE{ nday } = tmp.combinedData.UE;
    %barOn{ nday } = UE{ nday }.barOn - UE{ nday }.fixOn;
end

%% load un-aligned LFADS data and add in R to that Data
%loaddataPath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/data_masked_highCorr_rm/selected/';
loaddataPath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/highCorr_rm_noMasking/selected/';
%saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/data_masked_highCorr_rm/preAlignedData_added/selected/';
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/highCorr_rm_noMasking/preAlignedData_added/selected/';
if ~isdir(saveDir)
    mkdir(saveDir)
end

%% get aligned data
binsize_rescaled = 1;
for nday = 1:7
    barStart = tc_origin{nday}.UE.barOn - tc_origin{nday}.UE.fixOn;
    cueStart = tc_origin{nday}.UE.cueOn - tc_origin{nday}.UE.fixOn;
    targetStart = tc_origin{nday}.UE.targetOn - tc_origin{nday}.UE.fixOn;
    tc = [];
    for ntr = 1:numel(tc_origin{nday}.r.r)
        tc(ntr).spikes = tc_origin{nday}.r.r(ntr).spikes;
        tc( ntr ).barOnset = round( barStart( ntr ) /binsize_rescaled);
        tc( ntr ).cueOnset = round( cueStart( ntr ) / binsize_rescaled);
        tc( ntr ).targetStart = round( targetStart( ntr ) / binsize_rescaled);
    end

    % get PSTHs and rasters
    events = {'barOn', 'cueOn', 'targetOn'};
    conditions = struct;
    conditions.barOn = {'Vert', 'Hori'};
    conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
    window = round([-600  600] / binsize_rescaled);
    timePoints = window(1):window(2);

    %
    R = [];
    for nBar = 1:numel(conditions.barOn)
        trialIndicesPerCond.barOn.(conditions.barOn{nBar}) = tc_origin{nday}.UE.barType == nBar;
        % get R
        clear keepTrials_struct event_times tmp
        keepTrials_struct = tc(trialIndicesPerCond.barOn.(conditions.barOn{nBar}));
        event_times = [keepTrials_struct.barOnset];
        tmp(numel(keepTrials_struct)) = struct;
        for itrial = 1:numel(keepTrials_struct)
            tmp(itrial).spikeCounts = keepTrials_struct(itrial).spikes(:, event_times(itrial) + timePoints);
            tmp(itrial).condition = nBar;
            tmp(itrial).window = [-600, 600];
            tmp(itrial).type = ['barOn-' conditions.barOn{nBar}];
        end
        R = [R, tmp];
        
        for nCue = 1 : numel(conditions.cueOn)
            cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
            trialIndicesPerCond.cueOn.(cue_cond_field) = (tc_origin{nday}.UE.barType == nBar) & (tc_origin{nday}.UE.cueType == nCue);

            % get R
            clear keepTrials_struct event_times tmp
            keepTrials_struct = tc(trialIndicesPerCond.cueOn.(cue_cond_field));
            event_times = [keepTrials_struct.cueOnset];
            tmp(numel(keepTrials_struct)) = struct;
            for itrial = 1:numel(keepTrials_struct)
                tmp(itrial).spikeCounts = keepTrials_struct(itrial).spikes(:, event_times(itrial) + timePoints);
                tmp(itrial).condition = numel(conditions.barOn) + (nBar-1)*numel(conditions.cueOn) + nCue;
                tmp(itrial).window = [-600, 600];
                tmp(itrial).type = ['cueOn-' conditions.barOn{nBar} '-' conditions.cueOn{nCue}];
            end
            R = [R,tmp];
        end
    end
    mu_idx{nday} = 1:32;
    rm_channels{nday} = unique([hc_channels{nday}, drift_channels{nday}]);
    mu_idx{nday}(ismember(allChannels, rm_channels{nday})) = [];
    for i = 1:length(R)
        R(i).spikeCounts = R(i).spikeCounts(mu_idx{nday},:);
    end
    % add R to LFADS data
    clear combinedData
    loadName = [datasets(nday).shortName '_v2.mat'];
    fullDataName = fullfile(loaddataPath, loadName);
    load(fullDataName)
    combinedData.R = R;
    cd(saveDir);
    saveName = [datasets(nday).midName '_v3.mat'];
    save(saveName, 'combinedData', '-v7.3')
end

