%% define days and good neurons
dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
%dataset(5).date = '03272019'; % previously didn't work
dataset(5).date = '04062019';
dataset(6).date = '04252019';
dataset(7).date = '05022019';
dataset(8).date = '02082019';
dataset(9).date = '02132019';
dataset(10).date = '02142019';
dataset(11).date = '02152019';
dataset(12).date = '02162019';
dataset(13).date = '02262019';
dataset(14).date = '02282019';
dataset(15).date = '03012019';
dataset(16).date = '03022019';
dataset(17).date = '03032019';
dataset(18).date = '03312019';
dataset(19).date = '04012019';

%% good channel info
good_channels{1} = [8, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32];
good_channels{2} = [22, 23, 24, 26, 27, 28, 30, 31];
good_channels{3} = [14, 24, 25, 26, 30, 31];
good_channels{4} = [8, 21, 24, 25, 26, 27, 31, 32];
good_channels{5} = [17, 18, 19, 21, 23];
good_channels{6} = [9, 10, 21, 30, 31, 32];
good_channels{7} = [17, 24];
good_channels{8} = [11, 14, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26];
good_channels{9} = [4, 5, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 32];
good_channels{10} = [6, 15, 17, 18, 19, 22, 23, 24, 28, 29, 30, 31, 32];
good_channels{11} = [13, 16, 18, 21, 23, 25, 27, 29, 32];
good_channels{12} = [9, 11, 12, 15, 16, 19, 21, 23, 24, 29];
good_channels{13} = [11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 26, 27, 28, 29];
good_channels{14} = [13, 14, 15, 16, 23, 24, 25, 30, 32];
good_channels{15} = [10, 20, 30, 31];
good_channels{16} = [21, 24, 27, 28, 30, 32];
good_channels{17} = [18, 25, 27, 29, 30, 31];
good_channels{18} = [15 16 21 22 23 29 32];
good_channels{19} = [6 7 13 14 15 19];

%selected_days = [1, 2, 4, 8, 9, 10, 11, 13, 14, 17]; % first 10 best sessions
selected_days = [18, 19];

highCorr_channels{1} = 24;
highCorr_channels{2} = 24;
highCorr_channels{8} = [15, 22];
highCorr_channels{9} = [4, 17, 20];
highCorr_channels{10} = 19;
highCorr_channels{13} = 16;
highCorr_channels{18} = 22;
highCorr_channels{19} = 14;
highCorr_channels{numel(good_channels)} = [];

%%
for i = 1:numel(good_channels)
    keep_channels{i} = good_channels{i}(~ismember(good_channels{i}, highCorr_channels{i}));
end

%% load UE (all channels)
%baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/';
baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/withExtInp_lowerThresh/';
UE_baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/';
for iday = 1:numel(selected_days)
    nday = selected_days(iday);
    UE_file = ['UE_' dataset(nday).date '.mat'];
    load_UE_path = fullfile(UE_baseDir, UE_file);
    clear tmp
    tmp = load(load_UE_path);
    tmp.UE.date = dataset(selected_days(iday)).date;
    UE{ iday } = tmp.UE;
end

%% load spikes (no smooth or rebin) and only keep the channels that are responsive and not highly correlated
orginal_tc = {};
alf = {};
for iday = 1:numel(selected_days)
    nday = selected_days(iday);
    fileName = ['LIP_spiking_r_', dataset(nday).date, '.mat'];
    loadFile = fullfile(baseDir, fileName);
    disp( sprintf( 'loading day %g / %g', nday, numel( dataset ) ) );
    original_tc{iday} = load(loadFile);
    binsize_rescaled = 1;
    barStart = UE{iday}.barOn - UE{iday}.fixOn;
    cueStart = UE{iday}.cueOn - UE{iday}.fixOn;
    targetStart = UE{iday}.targetOn - UE{iday}.fixOn;
    %tc = [];
    for itr = 1:numel(original_tc{iday}.r.r)
        alf{iday}(itr).spikes = original_tc{iday}.r.r(itr).spikes(keep_channels{nday}, :);
        original_tc{iday}.r.r(itr).spikes = alf{iday}(itr).spikes;
        alf{iday}(itr).barOnset = round( barStart( itr ) /binsize_rescaled);
        alf{iday}(itr).cueOnset = round( cueStart( itr ) / binsize_rescaled);
        alf{iday}(itr).targetStart = round( targetStart( itr ) / binsize_rescaled);
    end
end

%% get aligned data
events = {'barOn', 'cueOn', 'targetOn'};
conditions = struct;
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
window = round([-50  350] / binsize_rescaled);
timePoints = window(1):window(2);
nConds = 16;
nTimes = length(timePoints);
pre_aligned = {};
for iday = 1:numel(selected_days)
    R = [];
    for nBar = 1:numel(conditions.barOn)
        for nCue = 1 : numel(conditions.cueOn)
            cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
            trialIndicesPerCond.cueOn.(cue_cond_field) = (UE{iday}.barType == nBar) & (UE{iday}.cueType == nCue);
            % get R
            clear keepTrials_struct event_times tmp
            keepTrials_struct = alf{iday}(trialIndicesPerCond.cueOn.(cue_cond_field));
            event_times = [keepTrials_struct.cueOnset];
            tmp(numel(keepTrials_struct)) = struct;
            for itrial = 1:numel(keepTrials_struct)
                tmp(itrial).spikeCounts = keepTrials_struct(itrial).spikes(:, event_times(itrial) + timePoints);
                tmp(itrial).condition = (nBar-1)*numel(conditions.cueOn) + nCue;
                tmp(itrial).window = [-50, 350];
                tmp(itrial).type = [conditions.barOn{nBar} '-' conditions.cueOn{nCue}];
            end
            R = [R,tmp];
        end
    end
    pre_aligned{iday}.R = R;
end

%% combine pre_aligned data and trialized continuous data and UE
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/newSignalProcessing/withPreAligned_withExtInp/lowerThresh/';
if ~isdir(saveDir)
    mkdir(saveDir)
end

for iday = 1:numel(selected_days)
    clear combinedData
    combinedData.r = Datasets.PulvinarTools.pulvinarData(original_tc{iday}.r.r);
    %combinedData.r = Datasets.PulvinarTools.pulvinarData_noExtInp(original_tc{iday}.r.r);
    combinedData.R = pre_aligned{iday}.R;
    combinedData.UE = UE{iday};
    combinedData.dayID = selected_days(iday);
    combinedData.channel_info = keep_channels{selected_days(iday)};
    combinedData.date = dataset(selected_days(iday)).date;
    cd(saveDir);
    saveName = [dataset(combinedData.dayID).date '_v2.mat'];
    save(saveName, 'combinedData', '-v7.3')
end

%%%
%day_seq = [8 9 10 11 1 13 14 17 2 4];
%loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/newSignalProcessing/withPreAligned_noExtInp/';
%for iday = 1:numel(day_seq)
%    clear combinedData
%    fname = sprintf( '%s%s_v1.mat', loadpath, dataset( iday ).date );
%    load( fname );
%    combinedData.channel_info = keep_channels{day_seq(iday)};
%    combinedData.date = dataset( iday ).date;
%    cd(loadpath);
%    saveName = [dataset(iday).date '_v2.mat'];
%    save(saveName, 'combinedData', '-v7.3')
%end
