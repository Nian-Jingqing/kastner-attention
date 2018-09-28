%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

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
% Pulvinar.Dataset(dc, '170320_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170324_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170327_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170329_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170331_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_lfadslite_20180408', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180408;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
par4( 1 ).c_in_factors_dim = 40;
par4( 1 ).c_keep_ratio = 0.5;
rc2.addParams( par4( 1 ) );
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

return;


%% load data and put them in Rstruct


for nData = 1:length(dc.datasets)
    realData(nData) = dc.datasets(nData).loadData();
    for trial = 1:length(realData(nData).R)
        realData(nData).R(trial).day = nData;
    end
end

% % use 170127 to test the code
% r_real = R.Rstruct(realData(1).R);
% r_realCopy = r_real.copy();


r_real = [realData.R];
r_real = R.Rstruct(r_real);
r_realCopy = r_real.copy();



%%
%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
datasetPath = ['/snel/share/share/derived/kastner/data_processed/singleSession/M20170608_PUL_all-g2-g3-g4-evokedSpiking/' ...
    'preAligned/CueOnArrayOnTargetDim_HoldRel/datasets'];

%% Locate and specify the datasets
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'preAligned_CO_AO_TD_HoldRel_20170608';

% add individual datasets
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;




%% load data and put it in Rstruct
r_real = dc.datasets(end).loadData(); % get the original dataset (for all neurons)
r_real = R.Rstruct(r_real.R); % put the dataset into R struct class
r_realCopy = r_real.copy();

%% pick up good neurons
nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ];

for itrial = 1: numel(r_realCopy.r)
    r_realCopy.r(itrial).spikeCounts = r_realCopy.r(itrial).spikeCounts(nIndices,:);
    r_realCopy.r(itrial).rfloc = r_realCopy.r(itrial).rfloc(nIndices,:);
    r_realCopy.r(itrial).day = 1;
end



%% get session info
% nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
binSize = 1;

%% smooth the firing rate
sigma = 50;
cutOff = sigma;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
r_realCopy.postSmoothCutOff( 'spike_smoothed', 'smoothed_cutOff', cutOff);

%% selecting alignment type (task period), trial type, and the days you want to include in the analysis
alignType_1 = 1;
alignType_2 = 2;
isHold_1 = 1;
isHold_2 = 1;
% nCond = 4;
days_vector = 1:6;


%% for each selected day, pull out the data for the day, extract the two periods for comparison, and do condition average
dayStruct(length(days_vector)).p1_condAvg = 1;
dayStruct(length(days_vector)).p2_condAvg = 1;

for d = 1: length(days_vector)
    r_thisDay = r_realCopy.r([r_realCopy.r.day] == days_vector(d)); % pull the day out
    r_p1 = r_thisDay([r_thisDay.alignType] == alignType_1 & [r_thisDay.isHoldTrial] == isHold_1); %pull the 1st period out
    r_p2 = r_thisDay([r_thisDay.alignType] == alignType_2 & [r_thisDay.isHoldTrial] == isHold_2); %pull the 2nd period out
    dayStruct(d).p1_condAvg = condAvg_general( r_p1, 'smoothed_cutOff', 'cueLoc' ); % compute condition avg for 1st period. nNeurons x (nConds x nTimes)
    dayStruct(d).p2_condAvg = condAvg_general( r_p2, 'smoothed_cutOff', 'cueLoc' ); % compute condition avg for 2nd period. nNeurons x (nConds x nTimes)
end

%% Concatenate all the days together nNeuronsTotal x (nConds x nTimes)
p1_condAvg = cat(1, dayStruct.p1_condAvg);
p2_condAvg = cat(1, dayStruct.p2_condAvg);

%% scale smoothed rate to firing rate
p1_condAvg = p1_condAvg * (1000/binSize);
p2_condAvg = p2_condAvg * (1000/binSize);

%% select time window for analysis for both time peorid
alignedTime_p1 = nTimesRaw/2 - cutOff;
startTime_p1 = alignedTime_p1 + 50 + 1;
endTime_p1 = startTime_p1 + 100 - 1;

alignedTime_p2 = nTimesRaw/2 - cutOff;
startTime_p2 = alignedTime_p2 + 50 + 1;
endTime_p2 = startTime_p2 + 100 - 1;

nTimesPerCond = nTimesRaw - 2*cutOff;

%% chop data into selected windows
nCond = length(unique([r_realCopy.r.cueLoc]));

for cond = 1:nCond
    startTimeThisCond_p1 = (cond-1)*nTimesPerCond + startTime_p1;
    endTimeThisCond_p1 = (cond-1)*nTimesPerCond + endTime_p1;
    condAvgStruct(cond).chopped_p1 = p1_condAvg( :, startTimeThisCond_p1 : endTimeThisCond_p1 );
    
    startTimeThisCond_p2 = (cond-1)*nTimesPerCond + startTime_p2;
    endTimeThisCond_p2 = (cond-1)*nTimesPerCond + endTime_p2;
    condAvgStruct(cond).chopped_p2 = p2_condAvg( :, startTimeThisCond_p2 : endTimeThisCond_p2 );
end

p1_condAvg_chopped = [condAvgStruct.chopped_p1];
p2_condAvg_chopped = [condAvgStruct.chopped_p2];

%% Normalize condition avg firing rate for specific analysis window

mean_p1_condAvg = mean(p1_condAvg_chopped, 2);
centered_p1_condAvg = bsxfun(@minus, p1_condAvg_chopped, mean_p1_condAvg);
centered_norm_p1_condAvg = centered_p1_condAvg./std(centered_p1_condAvg, 0, 2);
norm_p1_condAvg = bsxfun(@plus, centered_norm_p1_condAvg, mean_p1_condAvg);

mean_p2_condAvg = mean(p2_condAvg_chopped, 2);
centered_p2_condAvg = bsxfun(@minus, p2_condAvg_chopped, mean_p2_condAvg);
centered_norm_p2_condAvg = centered_p2_condAvg./std(centered_p2_condAvg, 0, 2);
norm_p2_condAvg = bsxfun(@plus, centered_norm_p2_condAvg, mean_p2_condAvg);

% %% 
% norm_p1_condAvg = p1_condAvg_chopped./(std(p1_condAvg_chopped')');
% norm_p2_condAvg = p2_condAvg_chopped./(std(p2_condAvg_chopped')');
% %% get rid of hi-corr neurons
% norm_p1_condAvg(3,:) = [];
% norm_p2_condAvg(3,:) = [];
% norm_p1_condAvg(4,:) = [];
% norm_p2_condAvg(4,:) = [];

%% try using centered normalized data
norm_p1_condAvg = centered_norm_p1_condAvg;
norm_p2_condAvg = centered_norm_p2_condAvg;

%% reshape normalized condition-average firing rate in the analysis window for each task period
nNeurons_total = size(norm_p1_condAvg, 1);
nTimesInWindow_p1 = endTime_p1 - startTime_p1 + 1;
nTimesInWindow_p2 = endTime_p2 - startTime_p2 + 1;
response_matrix_p1 = permute(reshape(norm_p1_condAvg, [nNeurons_total nTimesInWindow_p1 nCond]), [3 1 2]);
response_matrix_p2 = permute(reshape(norm_p2_condAvg, [nNeurons_total nTimesInWindow_p2 nCond]), [3 1 2]);



%% use the Tensor class to analyze this data
savedirOne = ['/snel/share/share/derived/kastner/nonLFADS_analysis/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/correlationMatrix/JanToMar/cueOn_arrayOn/'];
clear set

f1 = figure

t_1 = Tensor.Tensor;
t_1.t = response_matrix_p1;
subplot(2,2,1)
cmatrices_1 = t_1.calcCovarianceMatrices();
% mean_cmatrices_1 = mean(cmatrices_1.cy(:));
% cmatrices_1.cy(cmatrices_1.cy > 100) = mean_cmatrices_1;

imagesc( cmatrices_1.cy );
colorbar
title('unsorted matrix for cue-on');

% get an ordering for the neuron correlation matrix
order = Utils.getCrossCorrOrdering( cmatrices_1.cy, 50 );
subplot(2,2,2)
imagesc( cmatrices_1.cy( order, order ) )
colorbar
title('sorted matrix for cue-on');

% 
subplot(2,2,3)
t_2 = Tensor.Tensor;
t_2.t = response_matrix_p2;
cmatrices_2 = t_2.calcCovarianceMatrices();
imagesc( cmatrices_2.cy );
colorbar
title('unsorted matrix for array-on');

subplot(2,2,4)
imagesc( cmatrices_2.cy( order, order ) );
colorbar
title('sorted matrix for array-on based on cue-on ordering');

suptitle('CueOn vs. ArrayOn')
set(f1, 'Position', [204 41 1350 925]);
cd(savedirOne)
print(f1,'CueOn_vs_ArrayOn', '-dpng');

