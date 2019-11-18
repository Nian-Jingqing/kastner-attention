% this script calculate the condition avg activity across all the sessions and perform PCA to visualize lowD activity

%% define dates and good channels
%dataset(1).date = '02182019';
%dataset(2).date = '03062019';
%dataset(3).date = '03112019';
%dataset(4).date = '03142019';
%%%dataset(5).date = '03272019'; % previously didn't work
%dataset(5).date = '04062019';
%dataset(6).date = '04252019';
%dataset(7).date = '05022019';
%dataset(8).date = '02082019';
%dataset(9).date = '02132019';
%dataset(10).date = '02142019';
%dataset(11).date = '02152019';
%dataset(12).date = '02162019';
%dataset(13).date = '02262019';
%dataset(14).date = '02282019';
%dataset(15).date = '03012019';
%dataset(16).date = '03022019';
%dataset(17).date = '03032019';

% add in new sessions
%dataset(18).date = '03162019';
%dataset(19).date = '03312019';
%dataset(20).date = '04012019';
%dataset(21).date = '04052019';
%dataset(22).date = '04262019';

%good_channels{1} = [8, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32];
%good_channels{2} = [22, 23, 24, 26, 27, 28, 30, 31];
%good_channels{3} = [14, 24, 25, 26, 30, 31];
%good_channels{4} = [8, 21, 24, 25, 26, 27, 31, 32];
%good_channels{5} = [17, 18, 19, 21, 23];
%good_channels{6} = [9, 10, 21, 30, 31, 32];
%good_channels{7} = [17, 24];
%good_channels{8} = [11, 14, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26];
%good_channels{9} = [4, 5, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 32];
%good_channels{10} = [6, 15, 17, 18, 19, 22, 23, 24, 28, 29, 30, 31, 32];
%good_channels{11} = [13, 16, 18, 21, 23, 25, 27, 29, 32];
%good_channels{12} = [9, 11, 12, 15, 16, 19, 21, 23, 24, 29];
%good_channels{13} = [11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 26, 27, 28, 29];
%good_channels{14} = [13, 14, 15, 16, 23, 24, 25, 30, 32];
%good_channels{15} = [10, 20, 30, 31];
%good_channels{16} = [21, 24, 27, 28, 30, 32];
%good_channels{17} = [18, 25, 27, 29, 30, 31];
%good_channels{18} = [4, 8, 14, 16, 19, 21, 22];
%good_channels{19} = [15, 16, 21, 22, 23, 29, 32];
%good_channels{20} = [6, 7, 13, 14, 15, 19];
%good_channels{21} = [10, 14, 15, 17, 28, 32];
%good_channels{22} = [4, 21, 22, 23, 30, 31];

% lower s.d. more channels, didn't help
%good_channels2{1} = [2, 5, 8, 9, 11, 13, 18, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32]; %
%good_channels2{2} = [22, 23, 24, 26, 27, 28, 30, 31];
%good_channels2{3} = [9, 11, 14, 16, 24, 25, 27, 31, 32]; %
%good_channels2{4} = [8, 21, 24, 25, 26, 27, 31, 32];
%good_channels2{5} = [2, 14, 16, 19, 20, 23, 25, 30]; %
%good_channels2{6} = [9, 10, 21, 30, 31, 32];
%good_channels2{7} = [17, 24];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% below are the best sessions so far (11/18/2019) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset(1).date = '02082019';
dataset(2).date = '02132019';
dataset(3).date = '02142019';
dataset(4).date = '02152019';
dataset(5).date = '02182019';
dataset(6).date = '02262019';
dataset(7).date = '02282019';
dataset(8).date = '03032019';
dataset(9).date = '03062019';
dataset(10).date = '03142019';
%dataset(11).date = '03162019';
%dataset(12).date = '03312019';
%dataset(13).date = '04012019';
%dataset(14).date = '04052019';
%dataset(15).date = '04262019';

good_channels{1} = [11, 14, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26];
good_channels{2} = [4, 5, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 32];
good_channels{3} = [6, 15, 17, 18, 19, 22, 23, 24, 28, 29, 30, 31, 32];
good_channels{4} = [13, 16, 18, 21, 23, 25, 27, 29, 32];
good_channels{5} = [8, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32];
good_channels{6} = [11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 26, 27, 28, 29];
good_channels{7} = [13, 14, 15, 16, 23, 24, 25, 30, 32];
good_channels{8} = [18, 25, 27, 29, 30, 31];
good_channels{9} = [22, 23, 24, 26, 27, 28, 30, 31];
good_channels{10} = [8, 21, 24, 25, 26, 27, 31, 32];
%good_channels{11} = [4, 8, 14, 16, 19, 21, 22];
%good_channels{12} = [15, 16, 21, 22, 23, 29, 32];
%good_channels{13} = [6, 7, 13, 14, 15, 19];
%good_channels{14} = [10, 14, 15, 17, 28, 32];
%good_channels{15} = [4, 21, 22, 23, 30, 31];


% load UE and spiking activity, also smooth the spikes and rebin
% sigma is 100 and rebin at 10ms
% Directory of where the spiking data is stored
%spikingDataDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/';
spikingDataDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/withExtInp_lowerThresh/';
%spikingDataDir2 = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/lower_thresh/';
binsize_rescaled = 10;
[UE, alf] = loadUE_loadSpiking(spikingDataDir, dataset, good_channels, binsize_rescaled);
%[UE, alf2] = loadUE_loadSpiking(spikingDataDir2, dataset, good_channels2, binsize_rescaled);

%% calculate cond-avg matrix
events = {'barOn', 'cueOn', 'targetOn'};
conditions = struct;

% define all conditions to analyze
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};

% Define window to use to calculate cond-avg
window = round([-50, 350] / binsize_rescaled);
timePoints = window(1) : window(2);
nTimes = length(timePoints);
nConds = 16;

[cond_avg, channel_dayID] = getCondAvgAllDays(dataset, UE, alf, 'spikes', timePoints, nTimes, nConds, conditions, binsize_rescaled);
%[cond_avg2, channel_dayID2] = getCondAvgAllDays(dataset, UE, alf2, good_channels2, timePoints, nTimes, nConds, conditions, binsize_rescaled);

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
subplot(2,3,1)
plotCondAvgLowD(projected_all_days, 'All Days')
% plotCondAvgLowD is the function to plot a subplot with lowD activity of the conditions

% project each individual day
for nday = 1:5
    cond_avg_this_day = cond_avg(channel_dayID == nday, :);

    % regress each day's data to the global PCs to get a projection matrix
    %center the lowD activity (Not necessary actually, can removed)
    %dim_reduced_data_centered = bsxfun(@minus, dim_reduced_data, mean(dim_reduced_data, 2));
    % center the cond_avg for this day
    this_dataset_centered = bsxfun(@minus, cond_avg_this_day, mean(cond_avg_this_day, 2));
    % regress this day's data against global PCs to get a mapping
    %proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data_centered');
    proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data');
    % store the bias of this day, which is the which of each channel for this day
    proj_bias_this_day = mean(cond_avg_this_day, 2);

    % get projected individual day's lowD activity
    dim_reduced_data_this_day = proj_matrix_this_day' * this_dataset_centered;
    projected_this_day = permute(reshape(dim_reduced_data_this_day', [length(timePoints), 16, size(dim_reduced_data_this_day, 1)]), [1 3 2]);
    % reshape + permute to make dim_reduced_data a good shape for plotting,
    % which is nChannels x nTimes x nConditions

    % plot individual day
    subplot(2,3,nday+1)
    day_title = ['Day ' dataset(nday).date];
    plotCondAvgLowD(projected_this_day, day_title)
end
set(gcf, 'Position', [4 353 1914 733])
tightfig(gcf)

%% project individual days with more channels

for nday = [1 3 5]
    figure
    cond_avg_this_day = cond_avg2(channel_dayID2 == nday, :);

    % regress each day's data to the global PCs to get a projection matrix
    %center the lowD activity (Not necessary actually, can removed)
    %dim_reduced_data_centered = bsxfun(@minus, dim_reduced_data, mean(dim_reduced_data, 2));
    % center the cond_avg for this day
    this_dataset_centered = bsxfun(@minus, cond_avg_this_day, mean(cond_avg_this_day, 2));
    % regress this day's data against global PCs to get a mapping
    %proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data_centered');
    proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data');
    % store the bias of this day, which is the which of each channel for this day
    proj_bias_this_day = mean(cond_avg_this_day, 2);

    % get projected individual day's lowD activity
    dim_reduced_data_this_day = proj_matrix_this_day' * this_dataset_centered;
    projected_this_day = permute(reshape(dim_reduced_data_this_day', [length(timePoints), 16, size(dim_reduced_data_this_day, 1)]), [1 3 2]);
    % reshape + permute to make dim_reduced_data a good shape for plotting,
    % which is nChannels x nTimes x nConditions

    % plot individual day
    %subplot(2,4,nday+1)
    day_title = ['Day ' int2str(nday)];
    plotCondAvgLowD(projected_this_day, day_title)
end 