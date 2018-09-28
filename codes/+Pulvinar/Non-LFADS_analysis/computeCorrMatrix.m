%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')

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
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180214', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180214;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
rc2.addParams( par2 )
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

return;

%% load data and put them in Rstruct


for nData = 1:length(dc.datasets)
    realData(nData) = dc.datasets(nData).loadData();
end

% use 170127 to test the code
r_real = R.Rstruct(realData(1).R);
r_realCopy = r_real.copy();


% r_real = [realData.R];
% r_real = R.Rstruct(r_real);
% r_realCopy = r_real.copy();





% %% load data and put it in Rstruct
% r_real = dc.datasets(end).loadData(); % get the original dataset (for all neurons)
% r_real = R.Rstruct(r_real.R); % put the dataset into R struct class
% r_realCopy = r_real.copy();
% 
% %% pick up good neurons
% nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ];
% 
% for itrial = 1: numel(r_realCopy.r)
%     r_realCopy.r(itrial).spikeCounts = r_realCopy.r(itrial).spikeCounts(nIndices,:);
%     r_realCopy.r(itrial).rfloc = r_realCopy.r(itrial).rfloc(nIndices,:);
% end



%% get session info
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
binSize = 1;

%% smooth the firing rate
sigma = 50;
cutOff = sigma;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
r_realCopy.postSmoothCutOff( 'spike_smoothed', 'smoothed_cutOff', cutOff);

%% selecting alignment type and trial type, and define condition number, for each task peorid to analyze
alignType_1 = 1;
alignType_2 = 2;
isHold_1 = 1;
isHold_2 = 1;
nCond = 4;

%% calculate condition-averaged firing rate for each condition

% including: condition-average, allConditionAvgTogether, substract
% allAvgerage from each condition average
condAvgFiringRate_1 = condAvg_onSelectedAlignType( r_realCopy, alignType_1, isHold_1, binSize);
condAvgFiringRate_2 = condAvg_onSelectedAlignType( r_realCopy, alignType_2, isHold_2, binSize);



%% select time window for analysis for both time peorid
alignedTime_1 = nTimesRaw/2 - cutOff;
startTime_1 = alignedTime_1 + 40 + 1;
endTime_1 = startTime_1 + 300 - 1;

alignedTime_2 = nTimesRaw/2 - cutOff;
startTime_2 = alignedTime_2 + 40  + 1;
endTime_2 = startTime_2 + 300 - 1;

nTimesPerCond = nTimesRaw - 2*cutOff;

%% chop out the date for the selected time window and put them in to structure

% for task period 1
for cond = 1:nCond
    startTimeForCond_1 = (cond-1)*nTimesPerCond + startTime_1;
    endTimeForCond_1 = (cond-1)*nTimesPerCond + endTime_1;
    windowedCondStruct_1(cond).condAvgFiringRateInWindow = condAvgFiringRate_1( :, startTimeForCond_1 : endTimeForCond_1 );
end
condAvgFiringRateForAnalysis_1 = [windowedCondStruct_1.condAvgFiringRateInWindow];
% for task period 2
for cond = 1:nCond
    startTimeForCond_2 = (cond-1)*nTimesPerCond + startTime_2;
    endTimeForCond_2 = (cond-1)*nTimesPerCond + endTime_2;
    windowedCondStruct_2(cond).condAvgFiringRateInWindow = condAvgFiringRate_2( :, startTimeForCond_2 : endTimeForCond_2 );
end 
condAvgFiringRateForAnalysis_2 = [windowedCondStruct_2.condAvgFiringRateInWindow];


%% Normalize condition avg firing rate for specific analysis window
condAvgFiringRate_norm_1 = condAvgFiringRateForAnalysis_1./(std(condAvgFiringRateForAnalysis_1')');
condAvgFiringRate_norm_2 = condAvgFiringRateForAnalysis_2./(std(condAvgFiringRateForAnalysis_2')');

%% reshape normalized condition-average firing rate in the analysis window for each task period
nTimesInWindow_1 = endTime_1 - startTime_1 + 1;
nTimesInWindow_2 = endTime_2 - startTime_2 + 1;
response_matrix_1 = permute(reshape(condAvgFiringRate_norm_1, [nNeurons nTimesInWindow_1 nCond]), [3 1 2]);
response_matrix_2 = permute(reshape(condAvgFiringRate_norm_2, [nNeurons nTimesInWindow_2 nCond]), [3 1 2]);

%% use the Tensor class to analyze this data
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/CueOnArrayOnTargetDim_HoldRel/postAnalysis/' ...
    'withGoodNeurons_Run_20180201/Non-LFADS_analysis/correlationMatrix/withGoodNeurons/targetDim.vsTargetDim/'];

f1 = figure

t_1 = Tensor.Tensor;
t_1.t = response_matrix_1;
subplot(2,2,1)
cmatrices_1 = t_1.calcCovarianceMatrices();
% mean_cmatrices_1 = mean(cmatrices_1.cy(:));
% cmatrices_1.cy(cmatrices_1.cy > 100) = mean_cmatrices_1;

imagesc( cmatrices_1.cy );
colorbar
title('unsorted matrix for pre-targetDim');

% get an ordering for the neuron correlation matrix
order = Utils.getCrossCorrOrdering( cmatrices_1.cy, 20 );
subplot(2,2,2)
imagesc( cmatrices_1.cy( order, order ) )
colorbar
title('sorted matrix for pre-targetDim');

% 
subplot(2,2,3)
t_2 = Tensor.Tensor;
t_2.t = response_matrix_2;
cmatrices_2 = t_2.calcCovarianceMatrices();
imagesc( cmatrices_2.cy );
colorbar
title('unsorted matrix for post-targetDim');

subplot(2,2,4)
imagesc( cmatrices_2.cy( order, order ) );
colorbar
title('sorted matrix for post-targetDim based on pre-targetDim ordering');

suptitle('pre-TargetDim.vsPost-TargetDim')
set(f1, 'Position', [204 41 1350 925]);
cd(savedirOne)
% print(f1,'pre-TargetDim.vsPost-TargetDim', '-dpng');

