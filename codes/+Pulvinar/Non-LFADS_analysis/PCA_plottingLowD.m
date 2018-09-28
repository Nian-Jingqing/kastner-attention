%% set path to get tools for analysis and to get access to the dataset

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

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/runs'];
rc = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180201', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc.version = 20180201;

Pulvinar.definePulvinarRunParams;

return;

%% load data and put it in Rstruct
r_real = dc.datasets(end).loadData(); % get the original dataset (for all neurons)
r_real = R.Rstruct(r_real.R); % put the dataset into R struct class
r_realCopy = r_real.copy();
%% get session info
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
binSize = 1;

%% smooth the firing rate
sigma = 20;
cutOff = sigma;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
r_realCopy.postSmoothCutOff( 'spike_smoothed', 'smoothed_cutOff', cutOff);

%% selecting alignment type and trial type, and define condition number
alignType = 3;
isHold = 1;
nCond = 4;

%% calculate condition-averaged firing rate for each condition

% including: condition-average, allConditionAvgTogether, substract
% allAvgerage from each condition average
condAvgFiringRate = condAvg_onSelectedAlignType( r_realCopy, alignType, isHold, binSize);


%% Normalize condition avg firing rate 

%% select time window for analysis
alignedTime = nTimesRaw/2 - cutOff;
startTime = alignedTime + 280  + 1;
endTime = startTime + 300 - 1;
nTimesPerCond = nTimesRaw - 2*cutOff;


%% run PCA on selected analysis window
[coeff, score, latent, latent_tot, cum_variance] = PCA_onSelectedAnalysisWindow( condAvgFiringRate, startTime, endTime, nCond, nTimesPerCond );


%% plot condition average firing rate on the first couple of principle components

condColor = [[1 0 0]; [1 0.6 0]; [0 1 0]; [0 0 1]];
f1 = figure;

for condID = 1:4
    chopStart = (condID-1)*(endTime - startTime + 1) + 1;
    chopEnd = condID*(endTime - startTime + 1);
    
    scorePerCond = score( chopStart : chopEnd, : );
    
    sPlot_D1 = subplot(4,1,1)
    plot(scorePerCond(:, 1), 'color', condColor(condID,:));
    title(sPlot_D1, '1st component');
    hold on
    
    sPlot_D2 = subplot(4,1,2)
    plot(scorePerCond(:, 2), 'color', condColor(condID,:));
    title(sPlot_D2, '2nd component');
    hold on
    
    sPlot_D3 = subplot(4,1,3)
    plot(scorePerCond(:, 3), 'color', condColor(condID,:));
    title(sPlot_D3, '3rd component');
    hold on
    
    sPlot_D4 = subplot(4,1,4)
    plot(scorePerCond(:, 4), 'color', condColor(condID,:));
    title(sPlot_D4, '4th component');
    hold on
end
    
set(sPlot_D1, 'XTick',[1 200 400]);
set(sPlot_D1,'XTickLabels',{'-200','ArrayOnset','+200'}); 
set(sPlot_D2, 'XTick',[1 200 400]);
set(sPlot_D2,'XTickLabels',{'-200','ArrayOnset','+200'}); 
set(sPlot_D3, 'XTick',[1 200 400]);
set(sPlot_D3,'XTickLabels',{'-200','ArrayOnset','+200'}); 
set(sPlot_D4, 'XTick',[1 200 400]);
set(sPlot_D4,'XTickLabels',{'-200','ArrayOnset','+200'}); 

set(f1, 'Position', [680 97 779 989]);
    





%% multi period analysis

%% selecting alignment type and trial type, and define condition number
alignType1 = 1;
alignType2 = 2;
isHold = 1;
nCond = 4;

%% calculate condition-averaged firing rate for each condition

% including: condition-average, allConditionAvgTogether, substract
% allAvgerage from each condition average
condAvgFiringRate1 = condAvg_onSelectedAlignType( r_realCopy, alignType1, isHold, binSize);
condAvgFiringRate2 = condAvg_onSelectedAlignType( r_realCopy, alignType2, isHold, binSize);


%% select time window for analysis
alignedTime = nTimesRaw/2 - cutOff;
startTime1 = alignedTime + 40  + 1;
endTime1 = startTime1 + 300 - 1;

startTime2 = alignedTime + 40  + 1;
endTime2 = startTime2 + 300 - 1;

nTimesPerCond = nTimesRaw - 2*cutOff;


%% run PCA on selected analysis window
[coeff1, score1, latent1, latent_tot1, cum_variance1] = PCA_onSelectedAnalysisWindow( condAvgFiringRate1, startTime1, endTime1, nCond, nTimesPerCond );
[coeff2, score2, latent2, latent_tot2, cum_variance2] = PCA_onSelectedAnalysisWindow( condAvgFiringRate2, startTime2, endTime2, nCond, nTimesPerCond );

%% plot condition average firing rate on the first couple of principle components

condColor = [[1 0 0]; [1 0.6 0]; [0 1 0]; [0 0 1]];
f2 = figure;

for condID = 1:4
    chopStart = (condID-1)*nTimesPerCond + startTime2;
    chopEnd = (condID-1)*nTimesPerCond + endTime2;
    projectedCondAvgFiringRate = coeff1'*condAvgFiringRate2(:, chopStart:chopEnd);
    
    
    sPlot_D1 = subplot(4,1,1)
    plot(projectedCondAvgFiringRate(1, :), 'color', condColor(condID,:));
    title(sPlot_D1, '1st component');
    hold on
    
    sPlot_D2 = subplot(4,1,2)
    plot(projectedCondAvgFiringRate(2, :), 'color', condColor(condID,:));
    title(sPlot_D2, '2nd component');
    hold on
    
    sPlot_D3 = subplot(4,1,3)
    plot(projectedCondAvgFiringRate(3, :), 'color', condColor(condID,:));
    title(sPlot_D3, '3rd component');
    hold on
    
    sPlot_D4 = subplot(4,1,4)
    plot(projectedCondAvgFiringRate(4, :), 'color', condColor(condID,:));
    title(sPlot_D4, '4th component');
    hold on
end
    
% set(sPlot_D1, 'XTick',[1 200 400]);
% set(sPlot_D1,'XTickLabels',{'-200','ArrayOnset','+200'}); 
% set(sPlot_D2, 'XTick',[1 200 400]);
% set(sPlot_D2,'XTickLabels',{'-200','ArrayOnset','+200'}); 
% set(sPlot_D3, 'XTick',[1 200 400]);
% set(sPlot_D3,'XTickLabels',{'-200','ArrayOnset','+200'}); 
% set(sPlot_D4, 'XTick',[1 200 400]);
% set(sPlot_D4,'XTickLabels',{'-200','ArrayOnset','+200'}); 

set(f2, 'Position', [680 97 779 989]);

% ####################################### below is Chethan's demo
%% generate the data
% number of time samples
T = 100000; 
% low dimension
lowD = 2;
% high dimension
highD = 3;
% generate the lowD data using rand
X = randn(lowD, T);
%%
% make a projection matrix to highD space
W_X2Y = randn( lowD, highD);
% W_X2Y = [ 0 0 10; 0 5 0];

% project to highD space
Y = W_X2Y' * X;


%% actually run PCA
[coeff, score, latent] = pca( Y' );

% pca is finding projection matrix from high-D space to low-D space

% coeff is your projection matrix
% score is the representation of your data (Y) in the PC space
% latent is the amount of variance explained by each PC
%   must be normalized

% see how close the projected and real are
a = abs(score' - coeff'*Y);

%% 
% normalize the latents to get variance explained by each PC
latent_tot = latent / sum(latent);
% calculate the cumulative variance explained
cum_variance = cumsum( latent_tot );
