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
% add more datasets here if needed, same codec

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/runs'];
rc = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180125', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc.version = 20180125;

Pulvinar.definePulvinarRunParams;

return;
%% Post-running analysis - loading data and the output of LFADS
r_real = dc.datasets(end).loadData(); % get the original dataset (for all neurons)
r_real = R.Rstruct(r_real.R); % put the dataset into R struct class
for r_id = 1:length(rc.runs)
    run = rc.runs(r_id); % pull out run information
    run.loadSequenceData(); % load sequence data in that run
    run.loadPosteriorMeans(); % load posterior mean in that run
    run.addPosteriorMeansToSeq();
    r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct
end

%% set up nIndices to pull out the good neurons from r_real
r_realCopy = r_real.copy();
nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ];
neuronID = {'4a', '5b', '7b', '16b', '20a', '21b', '21c', '28a', '29a', '37a', '43a', '43c', '49b', '53a', '59a', '62a'};
% nIndices = [ 4 10 13 16 29 36 37 38 41 46 51 52 58 60 69 70 71 72 77 79 80 83 86 89 90 94 96 99 102 ];


%%
for itrial = 1: numel(r_realCopy.r)
    r_realCopy.r(itrial).spikeCounts = r_realCopy.r(itrial).spikeCounts(nIndices,:);
    r_realCopy.r(itrial).rfloc = r_realCopy.r(itrial).rfloc(nIndices,:);
end
    
    
%% get experiment info (nTrials, nTimes, nNeurons)
nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfads(1).r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfads(1).r(1).factors, 1);

%% Post-running analysis - Compare avg neuron firing rate and avg LFADS rate for different cue locations
% prepare for plotting
r_lfadsWhole = r_lfads(1).copy();
sigma = 50;
cutOff_smoothedSpiking = par.spikeBinMs*ceil(sigma/par.spikeBinMs);
r_realCopy.postSmoothCutOff('spikeCounts', 'rawSpike_cueOff', cutOff_smoothedSpiking);
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
r_realCopy.postSmoothCutOff( 'spike_smoothed', 'spike_cutOff', cutOff_smoothedSpiking);
cutOff_LFADSRates = ceil(sigma/par.spikeBinMs);
r_lfadsWhole.postSmoothCutOff( 'rates', 'rates_cutOff', cutOff_LFADSRates);
rbinned = r_realCopy.binData({'spike_cutOff'}, [par.spikeBinMs]);

AllTypeSmoothedSpikingMatrix = permute(cat(3, rbinned.spike_cutOff), [3 2 1]);
% re-organize the spiking data to nTrials x nTimes x nNeurons
AllTypeWholeLFADSRatesMatrix = permute(cat(3, r_lfadsWhole.r.rates_cutOff), [3 2 1]);
% re-organize the LFADS rates to nTrials x nTimes x nNeurons
AllTypeWholeSpikingMatrix = permute(cat(3, r_realCopy.r.rawSpike_cueOff), [3 2 1]);

nTimesLFADS = size(r_lfadsWhole.r(1).rates,2) - 2*cutOff_LFADSRates; 
% get new nTimesLFADS after cut off
nTimesRaw = size(r_realCopy.r(1).rawSpike_cueOff, 2);
AllclVector = arrayfun(@(x) x.cueLoc, r_realCopy.r); 
% pull out cue location for each trial and put into a vector
rfLoc = r_realCopy.r(1).rfloc;
% pull out receptive field for each neuron
%% separation by align types
alignType = 1; % select arrayOnset type
alignIx = [r_realCopy.r.alignType];
SmoothedSpikingMatrix = AllTypeSmoothedSpikingMatrix(alignIx == alignType,:,:);
WholeLFADSRatesMatrix = AllTypeWholeLFADSRatesMatrix(alignIx == alignType,:,:);
WholeSpikingMatrix = AllTypeWholeSpikingMatrix(alignIx == alignType,:,:);
clVector = AllclVector(alignIx == alignType);
nTrials = length(clVector);
%%

AllCueLoc = unique(clVector);
rebinSize = par.spikeBinMs; % data was rebinned from 1ms to 2ms during LFADS
AvgFiringRate = zeros(nNeurons, nTimesLFADS);
% initialize a avg firing rate matrix to store avg true firing rate (nNeurons x n Rebinned Times)
% for n = 1:nNeuron % loop over all neurons
AvgLFADSRate = zeros(nNeurons, nTimesLFADS);
% initialize a avg LFADS rate matrix to store avg LFADS rates (nNeurons x n Rebinned Times)
AvgFiringRate_InRF = zeros(nNeurons, nTimesLFADS);
AvgFiringRate_OffRF = zeros(nNeurons, nTimesLFADS);
AvgFiringRate_OutRF1 = zeros(nNeurons, nTimesLFADS);
AvgFiringRate_OutRF2 = zeros(nNeurons, nTimesLFADS);
% initialize avg Firing rates matrix for differet cue location
AvgLFADSRate_InRF = zeros(nNeurons, nTimesLFADS);
AvgLFADSRate_OffRF = zeros(nNeurons, nTimesLFADS);
AvgLFADSRate_OutRF1 = zeros(nNeurons, nTimesLFADS);
AvgLFADSRate_OutRF2 = zeros(nNeurons, nTimesLFADS);
% initialize avg LFADS rates matrix for differet cue location
TrialIndexInRF = false(nNeurons, nTrials);
TrialIndexOffRF = false(nNeurons, nTrials);
TrialIndexOutRF1 = false(nNeurons, nTrials);
TrialIndexOutRF2 = false(nNeurons, nTrials);
% initialize a matrix for each cue location to store the trial index. Ones
% would be index for trials that have the corresponding cue location

for n = 1:nNeurons % loop over all neurons

    AvgFiringRate(n,:) = (sum(SmoothedSpikingMatrix(:,:,n),1))*(1/nTrials)*(1000/par.spikeBinMs);
    % Store the avg firing rate to the avgFiringRate matrix
    AvgLFADSRate(n,:) = (sum(WholeLFADSRatesMatrix(:,:,n),1))*(1/nTrials);

    TrialIndexInRF(n,:) = clVector == rfLoc(n,1);
    nInRF = length(clVector(clVector == rfLoc(n,1)));
    % find the number of trials that cue loc is in RF
    TrialIndexOffRF(n,:) = clVector == rfLoc(n,2);
    % Find the index of the trials that cue loc is in Rf or off Rf
    nOffRF = length(clVector(clVector == rfLoc(n,2)));
    % find the number of trials that cue loc is Off RF
    OutRF = AllCueLoc(~ismember(AllCueLoc, rfLoc(n,:))); 
    % Find the cue locations that are not in or off Rf of this neuron
    TrialIndexOutRF1(n,:) = clVector == OutRF(1);
    % always find the index of the trials that cue loc number is smaller and
    % put these index into the OutRF1
    nOutRF1 = length(clVector(clVector == OutRF(1)));
    % find the number of trials that cue loc is OutRF 1
    TrialIndexOutRF2(n,:) = clVector == OutRF(2);
    % always find the index of the trials that cue loc number is larger and
    % put these index into the OutRF2
    nOutRF2 = length(clVector(clVector == OutRF(2)));
    % find the number of trials that cue loc is OutRF 2

    AvgFiringRate_InRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF)*(1000/par.spikeBinMs);
    AvgFiringRate_OffRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF)*(1000/par.spikeBinMs);
    AvgFiringRate_OutRF1(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1)*(1000/par.spikeBinMs);
    AvgFiringRate_OutRF2(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2)*(1000/par.spikeBinMs);

    AvgLFADSRate_InRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF);
    AvgLFADSRate_OffRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF);
    AvgLFADSRate_OutRF1(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1);
    AvgLFADSRate_OutRF2(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2);

end


%% plots of avg neuron firing rate and avg LFADS rate for different cue locations， plus spiking rasters
savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/postAnalysis/withGoodNeurons_Run_20180125/PSTH/Gaussian50ms/cueOnset/Ryan16/'];

cd(savedirOne);


for n = 1:nNeurons
    f1 = figure
    
    sPlot_avgRealRates = subplot(5,2,1) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
    plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
    hold on
    plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-748','CueOnset','+748'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgRealRates, 'Avg Real Firing Rate');
    ylabel('Firing Rate');
    set(gca,'xticklabel',{[]});
    
%     xlabel('time (ms)');
%     legend('show')
    
    sPlot_avgLFADSRates = subplot(5,2,2) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgLFADSRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgLFADSRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
    plot(AvgLFADSRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
    hold on
    plot(AvgLFADSRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-748','CueOnset','+748'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgLFADSRates, 'Avg LFADS Firing Rate');
    ylabel('Firing Rate');
%     set(gca,'YLim',sPlot_avgRealRates.YLim);
    xlabel('time (ms)');
    Legend = legend('show');
    set(Legend,'Position',[0.7387 0.6 0.1530 0.0655]);
    
    sPlot_InRFRaster = subplot(5,2,3)
    imagesc(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'});
    set(gca,'xticklabel',{[]});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    sPlot_OffRFRaster = subplot(5,2,5)
    imagesc(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'xticklabel',{[]})
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    % ylim([0 nOffRF]);
   
    
    sPlot_OutRF1Raster = subplot(5,2,7)
    imagesc(SmoothedSpikingMatrix(TrialIndexOutRF1(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'xticklabel',{[]})
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OutRF1Raster, 'Cue outside RF 1');
    ylabel('Trials');
    % ylim([0 nOutRF1]);
    
    
    sPlot_OutRF2Raster = subplot(5,2,9)
    imagesc(SmoothedSpikingMatrix(TrialIndexOutRF2(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-748','CueOnset','+748'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OutRF2Raster, 'Cue outside RF 2');
    ylabel('Trials');
    % ylim([0 nOutRF2]);
    xlabel('time (ms)');
    
    suptitle(['Neuron ' char(neuronID(n))]);
    set(f1, 'Position', [200 100 1000 1000]);
    
    print(f1,['Neuron ' char(neuronID(n))], '-dpng');
    close;
end

%% plots of avg neuron firing rate and avg LFADS rate for in and off RF cue locations， plus spiking rasters
savedirTwo = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7'...
    '/TargetDim/postAnalysis/AvgFiringPlusRasterOnOffRFCueLoc/Gaussian50ms'];

cd(savedirTwo);

n = 1;
for n = 1:20
    f2 = figure
    
    sPlot_avgRealRates = subplot(3,2,1) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    
   
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 2
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgRealRates, 'Avg Real Firing Rate');
    ylabel('Firing Rate');
    set(gca,'xticklabel',{[]});
    
    xlabel('time (ms)');
%    legend('show')
    
    sPlot_avgLFADSRates = subplot(3,2,2) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgLFADSRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgLFADSRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    
    
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgLFADSRates, 'Avg LFADS Firing Rate');
    ylabel('Firing Rate');
    set(gca,'YLim',sPlot_avgRealRates.YLim);
    xlabel('time (ms)');
    Legend = legend('show');
    set(Legend,'Position',[0.7387 0.5 0.1530 0.0655]);
    
    sPlot_InRFRaster = subplot(3,2,3)
    imagesc(WholeSpikingMatrix(TrialIndexInRF(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'});
    set(gca,'xticklabel',{[]});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    sPlot_OffRFRaster = subplot(3,2,5)
    imagesc(WholeSpikingMatrix(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesRaw nTimesRaw]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    % ylim([0 nOffRF]);
    xlabel('time (ms)');
    
    
    
    suptitle(['Neuron' int2str(n)]);
    set(f2, 'Position', [200 100 1000 900]);
    
    print(f2,['Neuron ' int2str(n)], '-dpng');
    close;
end
    
    

%% Plots of avg neuron firing rate and compare to avg LFADS rate （no separation between cueLoc）
savedirTwo = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/PSTH/PSTH_noCueLoc/Gaussian50ms'];

cd(savedirTwo);

for n = 1:20 % loop over all neurons
    f1 = figure
    subplot(2,1,1)
    plot(AvgLFADSRate(n,:), 'b', 'DisplayName','LFADS avg FR'); %Plot avg LFADS rate for that neuron
    hold on
    plot(AvgFiringRate(n,:), 'r', 'DisplayName','Real avg FR'); % plot avg true firing rate for that neuron
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-748','CueOnset','+748'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(['Neuron' int2str(n)]);
    ylabel('Firing Rate');
    xlabel('time (ms)');
    legend('show')
%     hold on
%     TopOfTimeZero = zeros(1,2);
%     TopOfTimeZero(1,1)= max(AvgLFADSRate(n,:)); TopOfTimeZero(1,2)= max(AvgFiringRate(n,:));
%     TopOfTimeZero = max(TopOfTimeZero);
%     plot([0.5*nTimesLFADS 0.5*nTimesLFADS],[1 TopOfTimeZero],'g', 'DisplayName','Time at TargetDim'); 
%     % Plot a vertical line at aligned time
    subplot(2,1,2)
    imagesc(SmoothedSpikingMatrix(:,:,n));
    set(gca,'XTick',[1 0.5*nTimesRaw nTimesRaw]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-748','CueOnset','+748'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title('Spiking raster');
    ylabel('Trials');
    % ylim([0 nOutRF2]);
    xlabel('time (ms)');
    
    suptitle(['Neuron' int2str(n)]);
    set(f1, 'Position', [200 100 1000 800]);
    
    
    print(f1,['Neuron ' int2str(n)], '-dpng');
    close;
end
    
    
