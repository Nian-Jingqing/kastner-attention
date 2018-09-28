%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools/')

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
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_PBTRun_20180403', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180403;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
%for nrun = 1:numel( par4 )
%    rc2.addParams( par4( nrun ) );
%end
par7.doPBT = true;
par7.PBTscript = '/snel/home/fzhu23/bin/PBT_HP_opt/pbt_opt/pbt_script_run_manager.py';
rc2.addParams( par7 )
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));





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

%% Post-running analysis - loading data and the output of LFADS
for nData = 1:length(dc.datasets)
    realData = dc.datasets(nData).loadData();
    r_real(nData) = R.Rstruct(realData.R);
end
% r_real = dc.datasets(1).loadData(); % get the original dataset (for all neurons)
% r_real = R.Rstruct(r_real.R); % put the dataset into R struct class

%% Run params 
% run_1:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.2;   c_kl_co_weight = 0.2
% run_2:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_3:    c_l2_gen_scale = 1;   c_kl_ic_weight = 0.8;   c_kl_co_weight = 0.8
% run_4:    c_l2_gen_scale = 10;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_5:    c_l2_gen_scale = 50;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5
% run_6:    c_l2_gen_scale = 100;   c_kl_ic_weight = 0.5;   c_kl_co_weight = 0.5



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
rbinned = r_realCopy.binData({'spikeCounts'}, [par7.spikeBinMs]);

%% load lfp data

myFolder = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/170127/';
fileName = '170127_cueOnArrayOnTargetDim_HoldRel_lfp.mat';
fullFileName = fullfile(myFolder, fileName);
data = load(fullFileName);
R = data.R;

%% filter the LFP data
filtHighCutoff = 15;
filtLowCutoff = 2;
Fs = 1000;

R = bandpassFilter( R, 'lfp', 'lfp_theta', filtHighCutoff, filtLowCutoff, Fs );




%% plotting the input with LFADS rate w/o real binned spiking
savedirBase = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withGoodNeurons_PBTRun_20180403/lfp/onlyRaster/'];
condNames = {'CueOnsetCueLoc1', 'CueOnsetCueLoc3', 'ArrayOnsetCueLoc1Hold', 'ArrayOnsetCueLoc3Hold', 'ArrayOnsetCueLoc1Rel', 'ArrayOnsetCueLoc3Rel', 'TargetDimCueLoc1', 'TargetDimCueLoc3' };
condIx = [r_lfadsWhole.r.conditionId];
chopStart = 1*(nTimesLFADS/4) + 1;
chopEnd = nTimesLFADS - 1*(nTimesLFADS/4);
chopStart_lfp = 301;
chopEnd_lfp = 1100;
channelLabel = cell(1, size(R(1).lfp, 1));
channelVector = size(R(1).lfp, 1) : -1 : 1;
for c = 1:size(R(1).lfp, 1)
    channelLabel{c} = int2str(channelVector(c));
end

div = 0.5;
clear set

for condType = 3:nCond
    trialsForThisCond = r_lfadsWhole.r(condIx == condType);
    trialsOfRealForThisCond = rbinned(condIx == condType);
    nTrialsThisCond = length(trialsForThisCond);
    R_selected = R(condIx == condType);
    if nTrialsThisCond > 40
        trialNum = 40;
    else
        trialNum = nTrialsThisCond;
    end
    for t = 1:trialNum
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


        % plot lfp data
%         lfpThisTrial = R_selected(t).lfp_theta;
%         chopped_lfp = lfpThisTrial(:,chopStart_lfp:chopEnd_lfp)';
%         dis = 0.5*(max(chopped_lfp(:)) - min(chopped_lfp(:)));
%         base = size(chopped_lfp, 2):-1:1;
%         riseBy = base*dis;
%         scrolled_lfp = chopped_lfp + riseBy;
%         sp(1) = subplot(2, 1, 1);
%         plot(scrolled_lfp, 'b', 'LineWidth', 1);
%         set(gca,'YTick', fliplr(riseBy));
%         set(gca,'YTickLabels',channelLabel);
%         set(gca,'XTick',[1 div*0.5*nTimesLFADS*par7.spikeBinMs div*nTimesLFADS*par7.spikeBinMs]);
%         set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
%         title(sp(1), 'LFP (2 - 15 Hz)')
%         ylabel('channels')
        
        % plot LFADS rate
        sp(2) = subplot(2, 1, 1);
        imagesc( lfadsRatesThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(2), 'LFADS rates');
        ylabel('Multi-units');
%         set(sp(7), 'FontSize', 7);
%         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
%         sp(8) = subplot(7,2,14);
%         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);


        % plot raw spiking
        sp(3) = subplot(2, 1, 2);
        imagesc( realSpikingThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(3), 'Real Spiking');
        ylabel('Multi-units');
% %         set(sp(8), 'FontSize', 7);
        
        
%         set(sp(1), 'Position', [0.1300 0.5456 0.7750 0.6])
%         set(sp(2), 'Position', [0.1300 0.1539 0.7750 0.3119]);
        suptitle(condNames{condType});
%         set(f1, 'Position', [279 53 648 913]);
%         set(f1, 'Position', [375 67 1079 899]);
        savedirOne = fullfile(savedirBase, condNames{condType});
        cd(savedirOne)
        print(f1, ['Trial ' int2str(t)], '-dpng');
        close;
    end
end
        
  
%%
        
          
            
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




