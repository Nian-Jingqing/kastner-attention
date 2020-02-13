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

[cond_avg_global] = getCondAvgAllDays_factors(dataset, UE, alf, timePoints, nTimes, nConds, conditions, binsize_rescaled);
%[cond_avg_eachDay, channel_dayID] = getCondAvgAllDays(dataset, UE, alf, 'spikes', timePoints, nTimes, nConds, conditions, binsize_rescaled); % for using smoothed spikes or rates
[cond_avg_eachDay, channel_dayID] = getCondAvgAllDays(dataset, UE, alf, 'factors', timePoints, nTimes, nConds, conditions, binsize_rescaled); % for using factors
%[cond_avg2, channel_dayID2] = getCondAvgAllDays(dataset, UE, alf2, good_channels2, timePoints, nTimes, nConds, conditions, binsize_rescaled);

%% do pca and get global PCs
all_data_means = mean( cond_avg_global, 2 ); % for using factors
%all_data_means = mean( cond_avg_eachDay, 2 ); % for using smoothed spikes or rates
all_data_centered = bsxfun(@minus, cond_avg_global, all_data_means); % for using factors
%all_data_centered = bsxfun(@minus, cond_avg_eachDay, all_data_means); % for using smoothed spikes or rates
[pca_proj_mat, pc_data, latent, tsquared, explained] = pca( all_data_centered', 'NumComponents', 10); % FZ changed from 10 to 40 on 12102019

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
    cond_avg_this_day = cond_avg_eachDay(channel_dayID == nday, :);

    % regress each day's data to the global PCs to get a projection matrix
    %center the lowD activity (Not necessary actually, can removed)
    %dim_reduced_data_centered = bsxfun(@minus, dim_reduced_data, mean(dim_reduced_data, 2));
    % center the cond_avg for this day
    % store the bias of this day, which is the which of each channel for this day
    proj_bias_this_day = mean(cond_avg_this_day, 2);
    this_dataset_centered = bsxfun(@minus, cond_avg_this_day, proj_bias_this_day);
    % regress this day's data against global PCs to get a mapping
    %proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data_centered');
    proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data');



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

%% project individual day's single-trial activity
figure
subplot(2,3,1)
plotCondAvgLowD(projected_all_days, 'All Days')
axis(gca, [-1 2 -2 2 -2 1])
for nday = 1:5
    cond_avg_this_day = cond_avg_eachDay(channel_dayID == nday, :);

    % regress each day's data to the global PCs to get a projection matrix
    %center the lowD activity (Not necessary actually, can removed)
    %dim_reduced_data_centered = bsxfun(@minus, dim_reduced_data, mean(dim_reduced_data, 2));
    % center the cond_avg for this day
    % store the bias of this day, which is the which of each channel for this day
    proj_bias_this_day = mean(cond_avg_this_day, 2);
    this_dataset_centered = bsxfun(@minus, cond_avg_this_day, proj_bias_this_day);
    % regress this day's data against global PCs to get a mapping
    %proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data_centered');
    proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data');


    % plot individual day
    subplot(2,3,nday+1)
    day_title = ['Day ' dataset(nday).date];
    %plotSingleTrialLowD(alf, UE, nday, proj_matrix_this_day, 'factors', timePoints, dataset(nday).date, proj_bias_this_day) % for using factors
    plotSingleTrialLowD(alf, UE, nday, proj_matrix_this_day, 'factors_smoothed', timePoints, dataset(nday).date, proj_bias_this_day, binsize_rescaled) % for using smoothed factors
    %plotSingleTrialLowD(alf, UE, nday, proj_matrix_this_day, 'rates', timePoints, dataset(nday).date, proj_bias_this_day, binsize_rescaled) % for using smoothed spikes
    %plotSingleTrialLowD_1(alf, UE, nday, pca_proj_mat, 'factors', timePoints, dataset(nday).date, all_data_means)
    axis(gca, [-1 2 -2 2 -2 1])
end
set(gcf, 'Position', [4 353 1914 733])
tightfig(gcf)

%% plot inidividual PCs vs time for all trials for all days
savedir = ['/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/' ...
           'postAnalysis/tc_newSP_PBT_191120/lowD_analysis/single_trials_meaningful/tmp_PCs_alltrialEachDay/'];
if ~isdir(savedir)
    mkdir(savedir);
end

cd(savedir)
cond{1} = 'TL-V';
cond{2} = 'BL-V';
cond{3} = 'TR-V';
cond{4} = 'BR-V';

pcID = 2;
for nday = 1:5
    figure
    for cueID = 1:4
        subplot(4,1,cueID)
        mean_pc = projected_all_days(:, pcID, cueID)';

        % perform this for each day
        cond_avg_this_day = cond_avg_eachDay(channel_dayID == nday, :);

        % regress each day's data to the global PCs to get a projection matrix
        %center the lowD activity (Not necessary actually, can removed)
        %dim_reduced_data_centered = bsxfun(@minus, dim_reduced_data, mean(dim_reduced_data, 2));
        % center the cond_avg for this day
        % store the bias of this day, which is the which of each channel for this day
        proj_bias_this_day = mean(cond_avg_this_day, 2);
        this_dataset_centered = bsxfun(@minus, cond_avg_this_day, proj_bias_this_day);
        % regress this day's data against global PCs to get a mapping
        %proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data_centered');
        proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data');
        plotSingleTrialPC(cueID, pcID, alf, UE, nday, proj_matrix_this_day, 'factors_smoothed', timePoints, dataset(nday).date, proj_bias_this_day, binsize_rescaled) % for using smoothed factors
        hold all;

        plot(mean_pc, 'Color', 'k', 'lineWidth', 3);
        hold on
        %legend('Location', 'best')
        title(cond{cueID})
        axis tight
    end
    set(gcf, 'Position', [21 4 843 1082])
    suptitle(['PC' int2str(pcID) ' Day' dataset(nday).date])
    saveFileName = ['PC' int2str(pcID) '-day' dataset(nday).date '-4conds'];
    print(gcf, saveFileName, '-dpng')
end

%% subtracting each PC's mean and concatenate data across days to get residual_tensor (nPCs x nTimes x nTrial) for each condition
% get trial indices to keep for each condition
for barID = 1:2
    for cueID = 1:8
        clear keepThisTrial trialsToKeep
        condID = (barID - 1)*8 + cueID;
        
        %cueID = 4;
        %barID = 1;
        field_to_analysis = 'factors_smoothed';
        mean_pc{condID} = projected_all_days(:, :, condID)';

        for nday = 1 : numel( alf )
            keepThisTrial{nday} = (UE{nday}.barType == barID) & (UE{nday}.cueType == cueID);
            trialsToKeep{nday} = find(keepThisTrial{nday});
        end

        totalTrialsToKeep = sum( cellfun( @sum, keepThisTrial ) );
        nPCs = size(pca_proj_mat, 2);
        residual_tensor{condID} = zeros(nPCs, nTimes, totalTrialsToKeep);
        lowD_tensor{condID} = zeros(nPCs, nTimes, totalTrialsToKeep);
        trialWhichDay{condID} = zeros(1, totalTrialsToKeep);

        tot_itrial = 1;
        for nday = 1:numel(alf)
            spikes_tensor{condID}{nday} = zeros(size(alf{nday}(1).spikes, 1), nTimes, numel(trialsToKeep{nday})); % initialize a spike tensor to store smoothed spiking
            rates_tensor{condID}{nday} = zeros(size(alf{nday}(1).rates_smoothed, 1), nTimes, numel(trialsToKeep{nday})); % initialize a spike tensor to store smoothed spiking
            
            cond_avg_this_day = cond_avg_eachDay(channel_dayID == nday, :);

            % regress each day's data to the global PCs to get a projection matrix
            %center the lowD activity (Not necessary actually, can removed)
            %dim_reduced_data_centered = bsxfun(@minus, dim_reduced_data, mean(dim_reduced_data, 2));
            % center the cond_avg for this day
            % store the bias of this day, which is the which of each channel for this day
            proj_bias_this_day = mean(cond_avg_this_day, 2);
            this_dataset_centered = bsxfun(@minus, cond_avg_this_day, proj_bias_this_day);
            % regress this day's data against global PCs to get a mapping
            %proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data_centered');
            proj_matrix_this_day = (this_dataset_centered' \ dim_reduced_data');
            for itr = 1:numel(trialsToKeep{nday})
                ntr = trialsToKeep{nday}(itr);
                data_this_trial = alf{nday}(ntr).(field_to_analysis)(:, alf{nday}(ntr).cueOnset + timePoints);
                data_this_trial_centered = bsxfun(@minus, data_this_trial, proj_bias_this_day);
                lowD_this_trial = proj_matrix_this_day' * data_this_trial_centered; % this should be nPCs x nTimes
                lowD_tensor{condID}(:, :, tot_itrial) = lowD_this_trial;
                residual_tensor{condID}(:, :, tot_itrial) = lowD_this_trial - mean_pc{condID};
                trialWhichDay{condID}(tot_itrial) = nday;
                spikes_tensor{condID}{nday}(:, :, itr) = alf{nday}(ntr).spikes(:, alf{nday}(ntr).cueOnset + timePoints);
                rates_tensor{condID}{nday}(:, :, itr) = alf{nday}(ntr).rates_smoothed(:, alf{nday}(ntr).cueOnset + timePoints);
                tot_itrial = tot_itrial + 1;
            end
        end

        if tot_itrial-1 ~= totalTrialsToKeep
            disp('WRONG WRONG WRONG')
        end

        % concatenate along the trials to get nPCs x nTimes*nTrials

        tmp = reshape(permute(residual_tensor{condID}, [2 3 1]), [], nPCs);
        residual_mat{condID} = tmp';
    end
end

%%
var_residual = zeros(nPCs, numel(residual_mat));
for iCond = 1:numel(residual_mat)
    var_residual(:, iCond) = var(residual_mat{iCond}, 0, 2);
end

%% plot the residual variance
cond{1} = 'TL-V-Ex';
cond{2} = 'BL-V-Ex';
cond{3} = 'TR-V-Ex';
cond{4} = 'BR-V-Ex';
cond{5} = 'TL-V-En';
cond{6} = 'BL-V-En';
cond{7} = 'TR-V-En';
cond{8} = 'BR-V-En';
cond{9} = 'TL-H-Ex';
cond{10} = 'BL-H-Ex';
cond{11} = 'TR-H-Ex';
cond{12} = 'BR-H-Ex';
cond{13} = 'TL-H-En';
cond{14} = 'BL-H-En';
cond{15} = 'TR-H-En';
cond{16} = 'BR-H-En';
figure
for iPC = 1:10
    plot(var_residual(iPC,:), '-o', 'DisplayName', ['PC' int2str(iPC)])
    hold on
end
xTickPos = 1:16;
set(gca, 'XTick', [xTickPos]);
set(gca, 'XTickLabels', cond);
legend('show')
legend('location', 'eastoutside')
ylabel('Variance')



%%
figure
for i = 1:10
    subplot(5,2,i)
    plot(squeeze(lowD_tensor{4}(i, :, :)));
    hold on
    plot(mean_pc{4}(i, :), 'Color', 'k', 'lineWidth', 3);
    title(['PC' int2str(i)])
    axis tight
    ylim([-1.5 2])
end
%suptitle('Condition TL-V')
suptitle('Condition BR-V')
set(gcf, 'Position', [207 4 865 1082]);


%%
figure
for i = 1:10
    subplot(5,2,i)
    plot(squeeze(residual_tensor{1}(i, :, :)), 'o');
    title(['PC' int2str(i)])
    axis tight
    ylim([-2 2])
end
%suptitle('Condition TL-V')
suptitle('Residuals, Condition BR-V')
set(gcf, 'Position', [207 4 865 1082]);


%% perform PCA on the residuals

all_residuals = [residual_mat{:}];

all_residual_means = mean( all_residuals, 2 ); % for using factors
%all_data_means = mean( cond_avg_eachDay, 2 ); % for using smoothed spikes or rates
all_residual_centered = bsxfun(@minus, all_residuals, all_residual_means); % for using factors
%all_data_centered = bsxfun(@minus, cond_avg_eachDay, all_data_means); % for using smoothed spikes or rates
[residual_pca_proj_mat, residual_pc_data, residual_latent, residual_tsquared, residual_explained] = pca( all_residual_centered', 'NumComponents', 10);    
    
%% perform single trial PCA on delay period
singleTrial_delay_lowD(alf, UE, binsize_rescaled, dataset)

%% project every trial in each condition and average across time points to get 1 number for each trials
for condID = 1:16
    trialPos{condID}.origin = zeros(1, size(residual_tensor{condID}, 3));
    for itr = 1:size(residual_tensor{condID}, 3)
        residual_this_trial = squeeze(residual_tensor{condID}(:,:,itr));
        residual_this_trial_centered = bsxfun(@minus, residual_this_trial, all_residual_means);
        dim_reduced_residual = residual_pca_proj_mat' * residual_this_trial_centered;
        trialPos{condID}.origin(itr) = mean(dim_reduced_residual(1, :));
    end
    trialPos{condID}.sorted = sort(trialPos{condID}.origin, 'ascend');
end

%% plot the distribution of trial position
figure
for condID = 1:16
    subplot(4,4,condID)
    histogram(trialPos{condID})
    title(['Cond' cond{condID}])
end

%% plot factor lowD, color order the trials by residual position,ALL DAYS
cond{1} = 'TL-V-Ex';
cond{2} = 'BL-V-Ex';
cond{3} = 'TR-V-Ex';
cond{4} = 'BR-V-Ex';
cond{5} = 'TL-V-En';
cond{6} = 'BL-V-En';
cond{7} = 'TR-V-En';
cond{8} = 'BR-V-En';
cond{9} = 'TL-H-Ex';
cond{10} = 'BL-H-Ex';
cond{11} = 'TR-H-Ex';
cond{12} = 'BR-H-Ex';
cond{13} = 'TL-H-En';
cond{14} = 'BL-H-En';
cond{15} = 'TR-H-En';
cond{16} = 'BR-H-En';
saveDir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/tc_newSP_PBT_191120/lowD_analysis/single_trials_meaningful/residual_ordered_factorsLowD/allFactorPCs/';
if ~isdir(saveDir)
    mkdir(saveDir)
end
cd(saveDir)
condToPlot = [1:4, 9:12];
cmap = colormap('parula');
nColors = size(cmap, 1);
figure
for iPC = 1:40
    for i = 1:8
        condID = condToPlot(i);
        subplot(4,2,i)
        lowD_this_cond = squeeze(lowD_tensor{condID}(iPC, :, :)); % this is nTimes x nTrials lowD activity for iPC
        for itr = 1:size(lowD_this_cond, 2)
            order = find(trialPos{condID}.sorted == trialPos{condID}.origin(itr));
            c = ceil(nColors * (order(1) / size(lowD_this_cond, 2)));
            plot(lowD_this_cond(:, itr), 'Color', cmap(c, :));
            hold on
        end
        plot(mean_pc{condID}(iPC, :), 'Color', 'k', 'lineWidth', 3);
        title(['Cond ' cond{condID}])
        axis tight
        %ylim([-1.5 2])
    end
    suptitle(['PC ' int2str(iPC)])
    set(gcf, 'Position', [207 4 865 1082]);        
    print(gcf, ['PC ' int2str(iPC)], '-dpng')        
end


%% plot smoothed spiking or rates, color order the trials by residual position, EACH DAY
baseDir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/tc_newSP_PBT_191120/lowD_analysis/single_trials_meaningful/residual_ordered_smoothedSpiking_100/';
%baseDir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/tc_newSP_PBT_191120/lowD_analysis/single_trials_meaningful/residual_ordered_smoothedRates/';
condToPlot = [1:4, 9:12];
cmap = colormap('parula');
nColors = size(cmap, 1);
for nday = 1:numel(alf)
    saveDir = fullfile(baseDir, dataset(nday).date);
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
    cd(saveDir)
    figure
    for iCh = 1:size(alf{nday}(1).spikes, 1) % for spikes
    %for iCh = 1:size(alf{nday}(1).rates_smoothed, 1)
        for i = 1:8
            condID = condToPlot(i);
            subplot(4,2,i)
            trialID_selected = find(trialWhichDay{condID} == nday);
            for itr = 1:numel(trialID_selected)
                ntr = trialID_selected(itr);
                order = find(trialPos{condID}.sorted == trialPos{condID}.origin(ntr));
                c = ceil(nColors * (order(1) / length(trialPos{condID}.origin)));
                ss_this_trial = squeeze(spikes_tensor{condID}{nday}(iCh, :, itr)); % for spikes
                %ss_this_trial = squeeze(rates_tensor{condID}{nday}(iCh, :, itr));
                plot(ss_this_trial, 'Color', cmap(c, :));
                hold on
            end
            plot(mean(spikes_tensor{condID}{nday}(iCh, :, :), 3), 'Color', 'k', 'lineWidth', 3); % for spikes
            %plot(mean(rates_tensor{condID}{nday}(iCh, :, :), 3), 'Color', 'k', 'lineWidth', 3);
            title(['Cond ' cond{condID}])
            axis tight
            %ylim([-1.5 2])
        end
        suptitle(['Neuron ' int2str(alf{nday}(1).channel_info(iCh))])
        set(gcf, 'Position', [207 4 865 1082]);        
        print(gcf, ['Neuron ' int2str(alf{nday}(1).channel_info(iCh))], '-dpng')
    end
end


%% plot factor lowD, color order the trials by residual position, EACH DAY
baseDir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/tc_newSP_PBT_191120/lowD_analysis/single_trials_meaningful/residual_ordered_factorsLowD/';
condToPlot = [1:4, 9:12];
cmap = colormap('parula');
nColors = size(cmap, 1);
for nday = 1:numel(alf)
    saveDir = fullfile(baseDir, dataset(nday).date);
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
    cd(saveDir)
    figure
    for iPC = 1:10
        for i = 1:8
            condID = condToPlot(i);
            subplot(4,2,i)
            trialID_selected = find(trialWhichDay{condID} == nday);
            for itr = 1:numel(trialID_selected)
                ntr = trialID_selected(itr);
                order = find(trialPos{condID}.sorted == trialPos{condID}.origin(ntr));
                c = ceil(nColors * (order(1) / length(trialPos{condID}.origin)));
                lowD_this_trial = squeeze(lowD_tensor{condID}(iPC, :, ntr));
                plot(lowD_this_trial, 'Color', cmap(c, :));
                hold on
            end
            plot(mean_pc{condID}(iPC, :), 'Color', 'k', 'lineWidth', 3);
            title(['Cond ' cond{condID}])
            axis tight
            %ylim([-1.5 2])
        end
        suptitle(['PC ' int2str(iPC)])
        set(gcf, 'Position', [207 4 865 1082]);        
        print(gcf, ['PC ' int2str(iPC)], '-dpng')        
    end
end


