% this script calculate the condition avg activity across all the sessions and perform PCA to visualize lowD activity

%% define dates and good channels
dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
%dataset(5).date = '03272019'; % previously didn't work
dataset(5).date = '04062019';
dataset(6).date = '04252019';
dataset(7).date = '05022019';

good_channels{1} = [8, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32];
good_channels{2} = [22, 23, 24, 26, 27, 28, 30, 31];
good_channels{3} = [14, 24, 25, 26, 30, 31];
good_channels{4} = [8, 21, 24, 25, 26, 27, 31, 32];
good_channels{5} = [17, 18, 19, 21, 23];
good_channels{6} = [9, 10, 21, 30, 31, 32];
good_channels{7} = [17, 24];

% load UE and spiking activity, also smooth the spikes and rebin
% sigma is 100 and rebin at 10ms
loadUE_loadSpiking

%% calculate cond-avg matrix
events = {'barOn', 'cueOn', 'targetOn'};
conditions = struct;

% define all conditions to analyze
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};

% define window to use to calculate cond-avg
window = round([-50, 350] / binsize_rescaled);
timePoints = window(1) : window(2);
nTimes = length(timePoints);
nConds = 16;

% cond_avg is the condition avg matrix (nChannel x nConditions*nTimes).
% There are 2 bar types and 8 cue types, so combined 16 conditions at cueOnset
% 16 is hard coded -- need to change
cond_avg = zeros(length([good_channels{:}]), nConds*nTimes);
% channel_dayID stores the day ID for each channel in the cond_avg matrix
channel_dayID = zeros(1, length([good_channels{:}]));

past_channelNum = 0; % how many channels in total in the past days
this_cond_where = 1:nTimes; % this is the inds of where the current cond is supposed to stored in the cond_avg matrix

% loop over each day and each bar type and each cue type
for nday = 1:numel(dataset)
    for nBar = 1:numel(conditions.barOn)        
        for nCue = 1 : numel(conditions.cueOn)
            % first select the condition to calculate 
            cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
            trialIndicesPerCond.cueOn.(cue_cond_field) = (UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue);

            % compute cond_avg of that condition. Used previous code when calculating PSTH
            [psth, raster_tensor, stderr] = preparePSTH_Manoj(alf{nday}, 'spikes', trialIndicesPerCond.cueOn.(cue_cond_field), [alf{nday}.cueOnset], timePoints, binsize_rescaled);
            cond_avg((1:length(good_channels{nday})) + past_channelNum, this_cond_where) = psth;
            this_cond_where = this_cond_where + nTimes;
            % update this_cond_where to indicate the inds for the next condition
        end
    end
    % fill in channel_dayID for the channels of the current day
    channel_dayID((1:length(good_channels{nday})) + past_channelNum) = nday;
    past_channelNum = past_channelNum + length(good_channels{nday});
    this_cond_where = 1:nTimes;
end

%% do pca and get global PCs
all_data_means = mean( cond_avg, 2 );
all_data_centered = bsxfun(@minus, cond_avg, all_data_means);
[pca_proj_mat, pc_data, latent, tsquared, explained] = pca( all_data_centered', 'NumComponents', 10);

%% plot projected activity

% project all-days data
dim_reduced_data = pca_proj_mat' * all_data_centered;
% dim_reduced_data is nChannels x nConditions*nTimes
projected_all_days = permute(reshape(dim_reduced_data', [nTimes, nConds, size(dim_reduced_data, 1)]), [1 3 2]);
% reshape + permute to make dim_reduced_data a good shape for plotting,
% which is nChannels x nTimes x nConditions

%projected = permute(reshape(pc_data, [length(timePoints), 16, size(pc_data, 2)]), [1 3 2]);
figure
subplot(2,4,1)
plotCondAvgLowD(projected_all_days, 'All Days')
% plotCondAvgLowD is the function to plot a subplot with lowD activity of the conditions

% project each individual day
for nday = 1:7
    cond_avg_this_day = cond_avg(channel_dayID == nday, :);

    % regress each day's data to the global PCs to get a projection matrix
    %center the lowD activity (Not necessary actually, can removed)
    dim_reduced_data_centered = bsxfun(@minus, dim_reduced_data, mean(dim_reduced_data, 2));
    % center the cond_avg for this day
    this_dataset_centered = bsxfun(@minus, cond_avg_this_day, mean(cond_avg_this_day, 2));
    % regress this day's data against global PCs to get a mapping
    proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data_centered');
    % store the bias of this day, which is the which of each channel for this day
    proj_bias_this_day = mean(cond_avg_this_day, 2);

    % get projected individual day's lowD activity
    dim_reduced_data_this_day = proj_matrix_this_day' * this_dataset_centered;
    projected_this_day = permute(reshape(dim_reduced_data_this_day', [length(timePoints), 16, size(dim_reduced_data_this_day, 1)]), [1 3 2]);
    % reshape + permute to make dim_reduced_data a good shape for plotting,
    % which is nChannels x nTimes x nConditions

    % plot individual day
    subplot(2,4,nday+1)
    day_title = ['Day ' int2str(nday)];
    plotCondAvgLowD(projected_this_day, day_title)
end
set(gcf, 'Position', [4 353 1914 733])