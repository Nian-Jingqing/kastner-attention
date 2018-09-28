%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/' ...
    'preAligned/multi-day_CoAoTdHoldRel_JanToMar/withGoodNeurons'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRel_JanToMar';

% add individual datasets
Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToMar/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180218', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180218;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
for nrun = 1:numel( par4 )
    rc2.addParams( par4( nrun ) );
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
for r_id = 1
    run = rc2.runs(r_id); % pull out run information
    run.loadSequenceData(); % load sequence data in that run
    run.loadPosteriorMeans(); % load posterior mean in that run
    run.addPosteriorMeansToSeq();
%     for nData = 1:length(run.sequenceData)
    nData = 1;
    RunID(r_id).r_lfads(nData) = R.Rstruct(run.sequenceData{nData});
%     end
%     r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct
end

%% Run params 
% run_1:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.2;   c_kl_co_weight = 0.2
% run_2:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_3:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.8;   c_kl_co_weight = 0.8
% run_4:    c_l2_gen_scale = 10;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_5:    c_l2_gen_scale = 50;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_6:    c_l2_gen_scale = 100;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5



%% Select the run and day you want to analyse
r_realCopy = r_real(6).copy();
r_lfadsWhole = RunID(1).r_lfads(6).copy();

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfadsWhole.r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfadsWhole.r(1).factors, 1);


%% pull out the Wrates


Day1_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170127_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day1_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170127_cueOnArrayOnTargetDim_HoldRel.h5_b:0' );

Day2_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170130_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day2_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170130_cueOnArrayOnTargetDim_HoldRel.h5_b:0' );

Day3_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170201_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day3_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170201_cueOnArrayOnTargetDim_HoldRel.h5_b:0' );

Day4_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170211_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day4_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170211_cueOnArrayOnTargetDim_HoldRel.h5_b:0' );

Day5_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170308_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day5_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170308_cueOnArrayOnTargetDim_HoldRel.h5_b:0' );

Day6_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170311_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day6_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/' ...
    'runs/withGoodNeurons_Run_20180218/param_gnmwq8/all/lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170311_cueOnArrayOnTargetDim_HoldRel.h5_b:0' );

%% Concatenate Wrates and b

Day1_Wrates = [Day1_Wrates, Day1_b];
Day2_Wrates = [Day2_Wrates, Day2_b];
Day3_Wrates = [Day3_Wrates, Day3_b];
Day4_Wrates = [Day4_Wrates, Day4_b];
Day5_Wrates = [Day5_Wrates, Day5_b];
Day6_Wrates = [Day6_Wrates, Day6_b];

%% Concatenate the days together

allDays_Wrates = [Day1_Wrates; Day2_Wrates; Day3_Wrates; Day4_Wrates; Day5_Wrates; Day6_Wrates];

%% set day identity
day1_id = ones(1, size(Day1_Wrates, 1));
day2_id = 2 * ones(1, size(Day2_Wrates, 1));
day3_id = 3 * ones(1, size(Day3_Wrates, 1));
day4_id = 4 * ones(1, size(Day4_Wrates, 1));
day5_id = 5 * ones(1, size(Day5_Wrates, 1));
day6_id = 6 * ones(1, size(Day6_Wrates, 1));

day_id = [day1_id, day2_id, day3_id, day4_id, day5_id, day6_id];


%% Normalize the Wrates matrix

norm_allDays_Wrates = allDays_Wrates./std(allDays_Wrates, 0, 2);


%% perform t-SNE

Wrates_tsne = tsne(norm_allDays_Wrates);

%% set up color vector
days_color = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0.8 0]; [0 1 1]; [1 0 1]];

%%
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
'multiDay_CO_AO_TD_HoldRel_JanToMar/postAnalysis/withGoodNeurons_Run_20180218/tSNE_onWrates'];
f1 = figure;
for nNeuron = 1:size(Wrates_tsne, 1)
    h = scatter(Wrates_tsne(nNeuron, 1), Wrates_tsne(nNeuron, 2));
    set(h, 'markerfacecolor', 'none');
%     set(h, 'markeredgecolor', days_color(day_id(nNeuron),:));
    set(h, 'markeredgecolor', 'b');
    hold on;
end
xlabel('1st dimension');
ylabel('2nd dimension');
    
cd(savedirOne);
print(f1, 'WithoutB', '-dpng');
close









