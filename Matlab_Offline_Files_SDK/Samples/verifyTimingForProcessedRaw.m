%% add dataset path
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')

%% set up day information
%datasets(1).shortName = '0208';
%datasets(1).longName = '02082019';
%datasets(2).shortName = '0218';
%datasets(2).longName = '02182019';
datasets(1).shortName = '0226';
datasets(1).longName = '02262019';
%datasets(4).shortName = '0227';
%datasets(4).longName = '02272019';
%datasets(5).shortName = '0308';
%datasets(5).longName = '03082019';
%datasets(6).shortName = '0310';
%datasets(6).longName = '03102019';
%datasets(7).shortName = '0311';
%datasets(7).longName = '03112019';

% assign rawDataRoot
rawDataRoot = '/snel/share/share/data/kastner/Manoj/PUL/';

%for day = 1:numel(datasets)
for day = 1

    %% construct the data path for loading
    folderName = ['Remy_' datasets(day).shortName '_PUL'];
    fileName = ['Remy_' datasets(day).longName '_PUL_MUA.mat'];
    fullFileName = fullfile(rawDataRoot, folderName, fileName);

    %% load the data
    S = load(fullFileName);

    % get the data into cell array. Each element would be the spike times for each MU

    % get rid of channels 1 - 4
    allChannels_spks = allChannels_spks(5:end);

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
end
