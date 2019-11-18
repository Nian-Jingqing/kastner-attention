% this script calculate the condition avg activity across all the sessions and perform PCA to visualize lowD activity

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

[cond_avg, channel_dayID] = getCondAvgAllDays(dataset, UE, alf, 'rates', timePoints, nTimes, nConds, conditions, binsize_rescaled);
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

%% project individual day's single-trial activity with more channels
figure
subplot(2,3,1)
plotCondAvgLowD(projected_all_days, 'All Days')
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

    % plot individual day
    subplot(2,3,nday+1)
    day_title = ['Day ' dataset(nday).date];
    plotSingleTrialLowD(alf, UE, nday, proj_matrix_this_day, 'rates', timePoints, dataset(nday).date)
end
set(gcf, 'Position', [4 353 1914 733])
tightfig(gcf)