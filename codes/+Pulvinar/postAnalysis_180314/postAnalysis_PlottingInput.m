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
for r_id = 2
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
%           c_co_dim = 6
% run_2:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.2;   c_kl_co_weight =
% 0.2
%           c_co_dim = 12


%% Select the run and day you want to analyse
r_realCopy = r_real(1).copy();
r_lfadsWhole = RunID(2).r_lfads(1).copy();

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfadsWhole.r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfadsWhole.r(1).factors, 1);
nInputs = size(r_lfadsWhole.r(1).controller_outputs, 1);
nCond = length(unique([r_lfadsWhole.r.conditionId]));

% %% loop over all the conditions to get condition-avg for all conditions
% condAvgWithBounds = [];
% for iCond = 1:nCond
%     [ condAvgWithBounds(iCond).condAvg, condAvgWithBounds(iCond).condAvg_up, condAvgWithBounds(iCond).condAvg_low ] = compute_condAvg( r_lfadsWhole, 'controller_outputs', iCond )
% end
% 
% % condAvgWithBounds is length = nCond structure array. Each element is one
% % condition
% % each field in each element is a nInput x nTimes matrix
%% loop over all the conditions and inputs
condStruct = [];
for c = 1:nCond
    condIx = [r_lfadsWhole.r.conditionId];
    trialsForThisCond = r_lfadsWhole.r(condIx == c);
    nTrialsForThisCond = length(trialsForThisCond);
    condStruct(c).inputTensor = permute(cat( 3, trialsForThisCond.controller_outputs), [ 3 2 1 ]);
    % extract out the trials for the condition and re-arrange the matrix
    % to nTrial x nTimes x nInputs
end

%% loop over all the conditions and inputs to plot


%%
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'postAnalysis/withGoodNeurons_Run_20180314/PlottingInputs/run2/FullTrials/'];
f1 = figure;
clear set
chopStart = 0*(nTimesLFADS/4) + 1;
chopEnd = nTimesLFADS - 0*(nTimesLFADS/4);
div = 1
for i = 1:nInputs
    
    sp((i-1)*3 + 1) = subplot(12, 3, (i-1)*3 + 1);
    c1_Input = condStruct(1).inputTensor(:, chopStart:chopEnd,i);
    shadedErrorBar([], c1_Input, {@mean, @(x) std(x)./sqrt(size(c1_Input, 1)) }, 'lineProps', '-r');
    hold on
    c2_Input = condStruct(2).inputTensor(:,chopStart:chopEnd,i);
    shadedErrorBar([], c2_Input, {@mean, @(x) std(x)./sqrt(size(c2_Input, 1)) }, 'lineProps', '-b');
    set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
    set(gca,'XTickLabels',{'-800','CueOnset','+800'});
    xlim([0 div*nTimesLFADS]);
    title(sp((i-1)*3 + 1), ['\fontsize{7}CueOnset - Input ' int2str(i)]);
    set(sp((i-1)*3 + 1), 'FontSize', 7);
    
    
    sp((i-1)*3 + 2) = subplot(12, 3, (i-1)*3 + 2);
    c3_Input = condStruct(3).inputTensor(:, chopStart:chopEnd,i);
    shadedErrorBar([], c3_Input, {@mean, @(x) std(x)./sqrt(size(c3_Input, 1)) }, 'lineProps', '-r');
    hold on
    c4_Input = condStruct(4).inputTensor(:,chopStart:chopEnd,i);
    shadedErrorBar([], c4_Input, {@mean, @(x) std(x)./sqrt(size(c4_Input, 1)) }, 'lineProps', '-b');
    set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
    set(gca,'XTickLabels',{'-800','ArrayOnset','+800'});
    xlim([0 div*nTimesLFADS]);
    title(sp((i-1)*3 + 2), ['\fontsize{7}ArrayOnset - Input ' int2str(i)]);
    set(sp((i-1)*3 + 2), 'FontSize', 7);
    
    
    
    sp((i-1)*3 + 3) = subplot(12, 3, (i-1)*3 + 3);
    c5_Input = condStruct(7).inputTensor(:, chopStart:chopEnd,i);
    h1(i) = shadedErrorBar([], c5_Input, {@mean, @(x) std(x)./sqrt(size(c5_Input, 1)) }, 'lineProps', '-r');
    hold on
    c6_Input = condStruct(8).inputTensor(:,chopStart:chopEnd,i);
    h2(i) = shadedErrorBar([], c6_Input, {@mean, @(x) std(x)./sqrt(size(c6_Input, 1)) }, 'lineProps', '-b');
    set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
    set(gca,'XTickLabels',{'-800','TargetDim','+800'});
    xlim([0 div*nTimesLFADS]);
    title(sp((i-1)*3 + 3), ['\fontsize{7}TargetDim - Input ' int2str(i)]);
    set(sp((i-1)*3 + 3), 'FontSize', 7);
end
    
legend([h1(1).mainLine, h2(1).mainLine], 'CueLoc 1', 'CueLoc 3', 'Location', 'best')
set(f1, 'Position', [25 40 1855 926]);
cd(savedirOne)
print(f1,'170127', '-dpng');
    
    


