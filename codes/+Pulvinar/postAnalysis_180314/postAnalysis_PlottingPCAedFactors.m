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
for r_id = 1
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

%% define whatever you need here
nPC_requested = 5;


%% Run params 
% run_1:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.2;   c_kl_co_weight = 0.2
%           c_co_dim = 6
% run_2:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.2;   c_kl_co_weight =
% 0.2
%           c_co_dim = 12

%% concatenate the factors to one matrix

factorsMatrix = [factorsAllDays.factors];

%% substract mean
all_factor_means = nanmean(factorsMatrix, 2);
all_factor_centered = bsxfun(@minus, factorsMatrix, all_factor_means);

%% actually run PCA
[pca_proj_mat, pc_data] = pca(all_factor_centered', 'NumComponents', nPC_requested);

%% loading data and the output of LFADS
%for r_id = 1:length(rc2.runs)
% r_id = 1
% run = rc2.runs(r_id); % pull out run information
% run.loadSequenceData(); % load sequence data in that run
% run.loadPosteriorMeans(); % load posterior mean in that run
% run.addPosteriorMeansToSeq();
allDayStruct = [];
for nData = 1:length(run.sequenceData)
    thisDay = R.Rstruct(run.sequenceData{nData});
    r_realCopy = r_real(1).copy();
    r_lfadsWhole = thisDay.copy();
    
    
    
    nTrials = length(r_realCopy.r); % get trial number
    nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
    nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
    nTimesLFADS = size(r_lfadsWhole.r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
    % modify this line if nTimes for different trials or runs are different.
    nFactors = size(r_lfadsWhole.r(1).factors, 1);
    nInputs = size(r_lfadsWhole.r(1).controller_outputs, 1);
    nCond = length(unique([r_lfadsWhole.r.conditionId]));
    
    
    
    condStruct = [];
    for c = 1:nCond
        condIx = [r_lfadsWhole.r.conditionId];
        trialsForThisCond = r_lfadsWhole.r(condIx == c);
        nTrialsForThisCond = length(trialsForThisCond);
        
        factors_thisCond = [trialsForThisCond.factors];
        factors_thisCond_centered = bsxfun(@minus, factors_thisCond, all_factor_means);
        dim_reduced_thisCond = pca_proj_mat' * factors_thisCond_centered;
        
        
        condStruct(c).factorTensor = permute(reshape(dim_reduced_thisCond, [ nPC_requested nTimesLFADS nTrialsForThisCond ]), [3 2 1]);
        % extract out the trials for the condition and re-arrange the matrix
        % to nTrial x nTimes x nPCs
    end
    
    allDayStruct(nData).factorByCond = condStruct;

end
%     r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct

%% plotting top 4 PCs for each condition
d = 12;
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'postAnalysis/withGoodNeurons_Run_20180314/factorsProjectedToPCs/170407/'];
f1 = figure;
clear set
chopStart = 0*(nTimesLFADS/4) + 1;
chopEnd = nTimesLFADS - 0*(nTimesLFADS/4);
div = 1
for i = 1:nPC_requested - 1
    
    sp((i-1)*3 + 1) = subplot(4, 3, (i-1)*3 + 1);
    c1_projected = [];
%     for d = 1:length(dc.datasets)
%         c1_projected = [c1_projected; allDayStruct(d).factorByCond(1).factorTensor(:, chopStart:chopEnd, i)];
%     end
    c1_projected = allDayStruct(d).factorByCond(1).factorTensor(:, chopStart:chopEnd, i);
    h1(i) = shadedErrorBar([], c1_projected, {@mean, @(x) std(x)./sqrt(size(c1_projected, 1)) }, 'lineProps', '-r');
    hold on
    c2_projected = [];
%     for d = 1:length(dc.datasets)
%         c2_projected = [c2_projected; allDayStruct(d).factorByCond(2).factorTensor(:, chopStart:chopEnd, i)];
%     end
    c2_projected = allDayStruct(d).factorByCond(2).factorTensor(:, chopStart:chopEnd, i);
    h2(i) = shadedErrorBar([], c2_projected, {@mean, @(x) std(x)./sqrt(size(c2_projected, 1)) }, 'lineProps', '-b');
    set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
    set(gca,'XTickLabels',{'-800','CueOnset','+800'});
    xlim([0 div*nTimesLFADS]);
    title(sp((i-1)*3 + 1), ['CueOnset - PC ' int2str(i)]);
    
    
    
    sp((i-1)*3 + 2) = subplot(4, 3, (i-1)*3 + 2);
    c3_projected = [];
%     for d = 1:length(dc.datasets)
%         c3_projected = [c3_projected; allDayStruct(d).factorByCond(3).factorTensor(:, chopStart:chopEnd, i)];
%     end
    c3_projected = allDayStruct(d).factorByCond(3).factorTensor(:, chopStart:chopEnd, i);
    h3(i) = shadedErrorBar([], c3_projected, {@mean, @(x) std(x)./sqrt(size(c3_projected, 1)) }, 'lineProps', '-r');
    hold on
    c4_projected = [];
%     for d = 1:length(dc.datasets)
%         c4_projected = [c4_projected; allDayStruct(d).factorByCond(4).factorTensor(:, chopStart:chopEnd, i)];
%     end
    c4_projected = allDayStruct(d).factorByCond(4).factorTensor(:, chopStart:chopEnd, i);
    h4(i) = shadedErrorBar([], c4_projected, {@mean, @(x) std(x)./sqrt(size(c4_projected, 1)) }, 'lineProps', '-b');
    
    hold on
    c5_projected = [];
    c5_projected = allDayStruct(d).factorByCond(5).factorTensor(:, chopStart:chopEnd, i);
    h5(i) = shadedErrorBar([], c5_projected, {@mean, @(x) std(x)./sqrt(size(c5_projected, 1)) }, 'lineProps', '-m');
    
    hold on
    c6_projected = [];
    c6_projected = allDayStruct(d).factorByCond(6).factorTensor(:, chopStart:chopEnd, i);
    h6(i) = shadedErrorBar([], c6_projected, {@mean, @(x) std(x)./sqrt(size(c6_projected, 1)) }, 'lineProps', '-c');
    
    set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
    set(gca,'XTickLabels',{'-800','ArrayOnset','+800'});
    xlim([0 div*nTimesLFADS]);
    title(sp((i-1)*3 + 2), ['ArrayOnset - PC ' int2str(i)]);
    
    
    
    sp((i-1)*3 + 3) = subplot(4, 3, (i-1)*3 + 3);
    c7_projected = [];
%     for d = 1:length(dc.datasets)
%         c7_projected = [c7_projected; allDayStruct(d).factorByCond(7).factorTensor(:, chopStart:chopEnd, i)];
%     end
    c7_projected = allDayStruct(d).factorByCond(7).factorTensor(:, chopStart:chopEnd, i);
    h7(i) = shadedErrorBar([], c7_projected, {@mean, @(x) std(x)./sqrt(size(c7_projected, 1)) }, 'lineProps', '-r');
    hold on
    c8_projected = [];
%     for d = 1:length(dc.datasets)
%         c6_projected = [c6_projected; allDayStruct(d).factorByCond(8).factorTensor(:, chopStart:chopEnd, i)];
%     end
    c8_projected = allDayStruct(d).factorByCond(8).factorTensor(:, chopStart:chopEnd, i);
    h8(i) = shadedErrorBar([], c8_projected, {@mean, @(x) std(x)./sqrt(size(c8_projected, 1)) }, 'lineProps', '-b');
    set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
    set(gca,'XTickLabels',{'-800','TargetDim','+800'});
    xlim([0 div*nTimesLFADS]);
    title(sp((i-1)*3 + 3), ['TargetDim - PC ' int2str(i)]);
    
end
    
lArray = legend([h3(1).mainLine, h4(1).mainLine, h5(1).mainLine, h6(1).mainLine], 'CueLoc 1 - hold', 'CueLoc 3 - hold', 'CueLoc 1 - rel', 'CueLoc 3 - rel', 'Location', 'southwest');
lArray.Position = [0.4281 0.75 0.0721 0.0622];
% lArray = legend([h3(1).mainLine, h4(1).mainLine], 'CueLoc 1 - hold', 'CueLoc 3 - hold', 'Location', 'best');
lCue = legend([h1(1).mainLine, h2(1).mainLine], 'CueLoc 1', 'CueLoc 3', 'Location', 'best')
lTarget = legend([h7(1).mainLine, h8(1).mainLine], 'CueLoc 1', 'CueLoc 3', 'Location', 'best')

lArray.FontSize = 6;
lCue.FontSize = 6;
lTarget.FontSize = 6;

set(f1, 'Position', [25 40 1855 926]);
cd(savedirOne)
print(f1,'4PC_FactorsProjection', '-dpng');





