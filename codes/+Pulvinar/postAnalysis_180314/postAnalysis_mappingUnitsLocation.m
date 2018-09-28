%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')

%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/preAligned/' ...
    'multi-day_CoAoTdHoldRel_JanToApr/withGoodNeurons_HoldRelSepForAO'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRelSepForAO_JanToApr';

% add individual datasets
Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170320_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170324_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170327_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170329_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170331_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180314', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180314;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
for nrun = 1:numel( par6 )
    rc2.addParams( par6( nrun ) );
end
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

return;


%% Post-running analysis - loading data and the output of LFADS
for nData = 1:length(dc.datasets)
    realData = dc.datasets(nData).loadData();
    r_real(nData) = R.Rstruct(realData.R);
end
% r_real = dc.datasets(1).loadData(); % get the original dataset (for all neurons)
% r_real = R.Rstruct(r_real.R); % put the dataset into R struct class




%% loading data and the output of LFADS
%for r_id = 1:length(rc2.runs)
factorsAllDays = [];
for r_id = 2
    run = rc2.runs(r_id); % pull out run information
    run.loadSequenceData(); % load sequence data in that run
    run.loadPosteriorMeans(); % load posterior mean in that run
    run.addPosteriorMeansToSeq();
    for nData = 1:length(run.sequenceData)
        currentDay = R.Rstruct(run.sequenceData{nData});
        factorsAllDays(nData).factors = [currentDay.r.factors];
        clear currentDay
    end
%     r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct
end


%% concatenate the factors to one matrix
factorsMatrix = [factorsAllDays.factors];
%% substract mean
all_factor_means = nanmean(factorsMatrix, 2);
all_factor_centered = bsxfun(@minus, factorsMatrix, all_factor_means);

%% actually run PCA
[pca_proj_mat, pc_data, latent] = pca(all_factor_centered');

% %% calculate total variance explained
% % normalize the latents to get variance explained by each PC
% latent_tot = latent / sum(latent);
% % calculate the cumulative variance explained
% cum_variance = cumsum( latent_tot );
% %% plot cum_variance
% 
% f1 = figure;
% dimension = 1:40;
% scatter(dimension, cum_variance)

%% get mapping from PCAed matrix to the original factors

W_p = pca_proj_mat';
W_pToFactors = inv(W_p' * W_p) * W_p';

%% pull out the Wrates


Day1_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170127_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day2_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170130_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day3_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170201_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day4_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170211_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day5_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170308_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day6_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170311_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day7_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170320_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day8_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170324_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day9_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170327_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day10_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170329_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day11_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170331_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day12_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170407_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

%% Concatenate the days together

allDays_Wrates = [Day1_Wrates; Day2_Wrates; Day3_Wrates; Day4_Wrates; Day5_Wrates; Day6_Wrates; Day7_Wrates; Day8_Wrates; Day9_Wrates; Day10_Wrates; Day11_Wrates; Day12_Wrates];

%% 

indices(1).goodNeurons = [1 3 4 6 7 8 11 13 14 15 16 17 19 20 21 22 23 24 25 26 29 30 31 32];
% 170127
indices(2).goodNeurons = [1 2 5 6 7 8 9 12 14 15 16 17 19 20 21 23 24 26 29 30 31];
% 170130
indices(3).goodNeurons = [1 5 8 10 12 13 27 30 32];
% 170201
indices(4).goodNeurons = [3 4 5 7 11 12 15 16 17 18 23 27 28 30];
% 170211
indices(5).goodNeurons = [9 15 18 19 20 28 30 31 32 33 34 37 38 39 42 46 49 55 58 62 63 64];
% 170308
indices(6).goodNeurons = [1 2 3 8 9 11 13 18 20 22 32 35 38 42 43 45 48 50 54 55 56 59];
% 170311
indices(7).goodNeurons = [4 5 9 12 16 20 24 25 29 30 31 34 35 37 38 40 41 44 49 50 51 52 53 54 56 57 59 61];
% 170320
indices(8).goodNeurons = [1 4 7 8 11 12 13 14 15 16 18 21 25 27 28 29 30 33 34 36 39 42 43 44 45 47 48 49 50 51 52 57 58 60 61 62 63 64];
% 170324
indices(9).goodNeurons = [1 2 3 4 5 10 11 13 15 16 18 19 21 22 24 25 27 31 32 35 39 40 41 42 43 44 45 47 50 52 53 54 55 56 57 58 59 60 61 62 63 64];
% 170327
indices(10).goodNeurons = [1 2 3 5 6 7 8 12 14 16 17 18 20 22 25 26 27 28 29 31 32 34 38 39 40 41 47 48 49 50 51 52 53 54 55 56 60];
% 170329
indices(11).goodNeurons = [7 9 10 11 13 14 17 24 25 29 30 31 32 34 36 37 41 46 47 50 51 52 55 57 58 59 62 64];
% 170331
indices(12).goodNeurons = [16 19 20 23 24 32];
% 170407

%% select the neurons that are in the subdivision

indices(1).ventral = 2:9;
indices(1).cortex = 10:32;
% 170127
indices(2).ventral = 16:25;
indices(2).cortex = 26:32;
% 170130
indices(3).ventral = 26:32;
indices(3).cortex = [];
% 170201
indices(4).ventral = 24:32;
indices(4).cortex = [];
% 170211
indices(5).ventral = [17:22, 50:56];
indices(5).cortex = [23:32, 57:64]
% 170308
indices(6).ventral = [8:9, 49:61];
indices(6).cortex = [10:32, 62:64];
% 170311
indices(7).ventral = [20:24, 43:51];
indices(7).cortex = [25:32, 52:64];
% 170320
indices(8).ventral = [11:18, 48:55];
indices(8).cortex = [19:32, 56:64];
% 170324
indices(9).ventral = [12:17, 57:64];
indices(9).cortex = 18:32;
% 170327
indices(10).ventral = 18:27;
indices(10).cortex = 28:32;
% 170329
indices(11).ventral = [7:12, 43:51];
indices(11).cortex = [13:32, 52:64];
% 170331
indices(12).ventral = 7:12;
indices(12).cortex = 13:32;
% 170407


% indices(1).ventral = 2:9;
% indices(1).cortex = 28:32;
% % 170127
% indices(2).ventral = 16:25;
% indices(2).cortex = 28:32;
% % 170130
% indices(3).ventral = 26:32;
% indices(3).cortex = [];
% % 170201
% indices(4).ventral = 24:32;
% indices(4).cortex = [];
% % 170211
% indices(5).ventral = [17:22, 50:56];
% indices(5).cortex = [28:32, 60:64]
% % 170308
% indices(6).ventral = [8:9, 49:61];
% indices(6).cortex = [28:32, 62:64];
% % 170311
% indices(7).ventral = [20:24, 43:51];
% indices(7).cortex = [28:32, 60:64];
% % 170320
% indices(8).ventral = [11:18, 48:55];
% indices(8).cortex = [28:32, 60:64];
% % 170324
% indices(9).ventral = [12:17, 57:64];
% indices(9).cortex = 28:32;
% % 170327
% indices(10).ventral = 18:27;
% indices(10).cortex = 28:32;
% % 170329
% indices(11).ventral = [7:12, 43:51];
% indices(11).cortex = [28:32, 60:64];
% % 170331
% indices(12).ventral = 7:12;
% indices(12).cortex = 28:32;
% % 170407



indices(1).dorsal = 1;
% 170127
indices(2).dorsal = 6:15;
% 170130
indices(3).dorsal = 16:25;
% 170201
indices(4).dorsal = 14:23;
% 170211
indices(5).dorsal = [7:16, 40:49];
% 170308
indices(6).dorsal = [1:7, 39:48];
% 170311
indices(7).dorsal = [10:19, 33:42];
% 170320
indices(8).dorsal = [1:10, 38:47];
% 170324
indices(9).dorsal = [2:11, 47:56];
% 170327
indices(10).dorsal = 8:17;
% 170329
indices(11).dorsal = [1:6, 33:42];
% 170331
indices(12).dorsal = 1:6;
% 170407

%% index the neurons with labels of ventral, dorsal or cortex

for d = 1:12
    [~, ventralInGood] = intersect(indices(d).goodNeurons, indices(d).ventral);
    indices(d).goodIsVentral = false(size(indices(d).goodNeurons));
    indices(d).goodIsVentral(ventralInGood) = true;
    [~, dorsalInGood] = intersect(indices(d).goodNeurons, indices(d).dorsal);
    indices(d).goodIsDorsal = false(size(indices(d).goodNeurons));
    indices(d).goodIsDorsal(dorsalInGood) = true;
    [~, cortexInGood] = intersect(indices(d).goodNeurons, indices(d).cortex);
    indices(d).goodIsCortex = false(size(indices(d).goodNeurons));
    indices(d).goodIsCortex(cortexInGood) = true;
end

%% concatenate the indices from all days together
isVentral = [indices.goodIsVentral];
isDorsal = [indices.goodIsDorsal];
isCortex = [indices.goodIsCortex];


%% set day identity
day1_id = ones(1, size(Day1_Wrates, 1));
day2_id = 2 * ones(1, size(Day2_Wrates, 1));
day3_id = 3 * ones(1, size(Day3_Wrates, 1));
day4_id = 4 * ones(1, size(Day4_Wrates, 1));
day5_id = 5 * ones(1, size(Day5_Wrates, 1));
day6_id = 6 * ones(1, size(Day6_Wrates, 1));
day7_id = 7 * ones(1, size(Day7_Wrates, 1));
day8_id = 8 * ones(1, size(Day8_Wrates, 1));
day9_id = 9 * ones(1, size(Day9_Wrates, 1));
day10_id = 10 * ones(1, size(Day10_Wrates, 1));
day11_id = 11 * ones(1, size(Day11_Wrates, 1));
day12_id = 12 * ones(1, size(Day12_Wrates, 1));


day_id = [day1_id, day2_id, day3_id, day4_id, day5_id, day6_id, day7_id, day8_id, day9_id, day10_id, day11_id, day12_id];

%% calculate W_pToRate

W_pToRate = allDays_Wrates * W_pToFactors;

%% Select how many PCs you want for t-SNE analysis

nPCs = 4;


%% Normalize the W_pToRate matrix
mean_W_pToRate = mean(W_pToRate(:, 1:nPCs), 2);
centered_W_pToRate = bsxfun(@minus, W_pToRate(:, 1:nPCs), mean_W_pToRate);
centered_norm_W_pToRate = centered_W_pToRate./std(centered_W_pToRate, 0, 2);
norm_W_pToRate = bsxfun(@plus, centered_norm_W_pToRate, mean_W_pToRate);


%% perform t-SNE

W_tsne = tsne(norm_W_pToRate);




%% set up color vector
days_color = [[1 0 0]; [1 0.8 0]; [1 1 0 ]; [0.5 1 0]; [0 1 0]; [0 1 0.2]; [0 1 0.8]; [0 1 1]; [0 0 1]; [0.2 0 1]; [0.8 0 1]; [1 0 1]];


%% Plot clusters seperated by recording locations

savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/' ...
    'withGoodNeurons_Run_20180314/tSNE_onWratesWithFactorPCs/'];
f1 = figure;
clear set

h_ventral = scatter(W_tsne(isVentral, 1), W_tsne(isVentral, 2), 'DisplayName','ventral');
set(h_ventral, 'markerfacecolor', 'none');
set(h_ventral, 'markeredgecolor', 'b');
% set(h_ventral, 'Legend', 'b');
% legend(h_ventral, {'ventral'});
hold on;

h_dorsal = scatter(W_tsne(isDorsal, 1), W_tsne(isDorsal, 2), 'DisplayName','dorsal');
set(h_dorsal, 'markerfacecolor', 'none');
set(h_dorsal, 'markeredgecolor', 'r');
% set(h_dorsal, 'Legend', 'r');
% legend(h_dorsal, {'dorsal'});
hold on;

h_cortex = scatter(W_tsne(isCortex, 1), W_tsne(isCortex, 2), 'DisplayName','cortex');
set(h_cortex, 'markerfacecolor', 'none');
set(h_cortex, 'markeredgecolor', 'g');
% set(h_cortex, 'Legend', 'g');
% legend(h_cortex, {'cortex'});
hold on;

xlabel('1st dimension');
ylabel('2nd dimension');
legend
    
cd(savedirOne);
print(f1, '4PCs_codedByRecordingLocations', '-dpng');
% close





%%
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/' ...
    'withGoodNeurons_Run_20180314/tSNE_onWratesWithFactorPCs/'];
f1 = figure;
for nNeuron = 1:size(W_tsne, 1)
    h = scatter(W_tsne(nNeuron, 1), W_tsne(nNeuron, 2));
    set(h, 'markerfacecolor', 'none');
    set(h, 'markeredgecolor', days_color(day_id(nNeuron),:));
%     set(h, 'markeredgecolor', 'b');
    hold on;
end
xlabel('1st dimension');
ylabel('2nd dimension');
    
cd(savedirOne);
% print(f1, '10PCs_codedByDays', '-dpng');
% close

