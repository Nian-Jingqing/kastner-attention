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
%           c_co_dim = 6
% run_2:    c_l2_gen_scale = 100;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
%           c_co_dim = 6



%% Select the run and day you want to analyse
r_realCopy = r_real(1).copy();
r_lfadsWhole = RunID(1).r_lfads(1).copy();

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfadsWhole.r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfadsWhole.r(1).factors, 1);
nInputs = size(r_lfadsWhole.r(1).controller_outputs, 1);
nCond = length(unique([r_lfadsWhole.r.conditionId]));


%% 
rbinned = r_realCopy.binData({'spikeCounts'}, [par4.spikeBinMs]);


%% plotting the input with LFADS rate w/o real binned spiking
savedirBase = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToMar/postAnalysis/' ...
    'withGoodNeurons_Run_20180218/InputsWithLFADSRates/withRealFiring_400Prior400Post/'];
condNames = {'CueOnsetCueLoc1', 'CueOnsetCueLoc3', 'ArrayOnsetCueLoc1', 'ArrayOnsetCueLoc3', 'TargetDimCueLoc1', 'TargetDimCueLoc3' };
condIx = [r_lfadsWhole.r.conditionId];
chopStart = 1*(nTimesLFADS/4) + 1;
chopEnd = nTimesLFADS - 1*(nTimesLFADS/4);
div = 0.5
clear set
for condType = 1:nCond
    trialsForThisCond = r_lfadsWhole.r(condIx == condType);
    trialsOfRealForThisCond = rbinned(condIx == condType);
    nTrialsThisCond = length(trialsForThisCond);
    for t = 1:40
        lfadsRatesThisTrial = trialsForThisCond(t).rates;
        inputThisTrial = trialsForThisCond(t).controller_outputs;
        realSpikingThisTrial = trialsOfRealForThisCond(t).spikeCounts;
        f1 = figure;
%         for i = 1:nInputs
%             sp(i) = subplot(7, 2, i*2 - 1);
%             plot(inputThisTrial(i, chopStart:chopEnd), 'b');
%             set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
%             set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
%             set(sp(i), 'FontSize', 7);
%             xlim([0 div*nTimesLFADS]);
%             title(sp(i), ['Input ' int2str(i)]);
%         end
%         sp(7) = subplot(7,2,13);
%         sp(7) = subplot('position', [0.1300 0.0339 0.3347 0.150]);
        sp(1) = subplot(2, 1, 1)
        imagesc( lfadsRatesThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(1), 'LFADS rates');
        ylabel('units');
%         set(sp(7), 'FontSize', 7);
%         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
%         sp(8) = subplot(7,2,14);
%         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);
        sp(2) = subplot(2, 1, 2)
        imagesc( realSpikingThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(2), 'Real Spiking');
        ylabel('units');
%         set(sp(8), 'FontSize', 7);
        
        suptitle(condNames{condType});
%         set(f1, 'Position', [279 53 648 913]);
%         set(f1, 'Position', [400 38 1219 928]);
        savedirOne = fullfile(savedirBase, condNames{condType});
        cd(savedirOne)
        print(f1, ['Trial ' int2str(t)], '-dpng');
        close;
    end
end
        
  
%%
        
          
            
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




