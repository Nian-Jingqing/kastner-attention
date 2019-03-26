%% add your paths here.

% add paths for Feng
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/postAnalysisCodes')

%% test jPCA code
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools/fromMarksLibraries')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools/CircStat2010d')

%% build run
buildRuns_20181023

%% load LFADS rates, factors, real spiking, smoothed spiking, UEs into cell array for each day
% please specify smoothing window here
sigma = 10;
binsize = par.spikeBinMs;
rfactor = 1;
binsize_rescaled = binsize * rfactor;
loadAligned_twoLocations

%% select alignment type you want to plot
% cueOnset, arrayOnset, or dim?
alignType = 1; % cueAligned
%%

AllTypeSmoothedSpikingMatrix = permute(cat(3, alf{nday}.binned_spikes), [3 2 1]);
% % re-organize the spiking data to nTrials x nTimes x nNeurons
AllTypeWholeLFADSRatesMatrix = permute(cat(3, alf{nday}.rates), [3 2 1]);
% % re-organize the LFADS rates to nTrials x nTimes x nNeurons
% % AllTypeWholeSpikingMatrix = permute(cat(3, r_realCopy.r.rawSpike_cueOff), [3 2 1]);

% nTimesLFADS = size(r_lfadsWhole.r(1).rates,2) - 2*cutOff_LFADSRates; 
% % get new nTimesLFADS after cut off
% nTimesRaw = size(r_realCopy.r(1).rawSpike_cueOff, 2);

% AllclVector = arrayfun(@(x) x.cueLoc, r_realCopy.r); 
% pull out cue location for each trial and put into a vector
% rfLoc = r_realCopy.r(1).rfloc;
% pull out receptive field for each neuron
%% separation by align types
% alignType = 1; % select arrayOnset type
if alignType ~= 3
    alignInds = ( alignType - 1 ) * numel( UEs{ nday }.cueLoc ) + (1 : numel( UEs{ nday }.cueLoc ));
    clVector = UEs{ nday }.cueLoc;
else
    alignInds = ( alignType - 1 ) * numel( UEs{ nday }.cueLoc ) + (1 : numel( find(UEs{ nday }.isHoldTrial )));
    clVector = UEs{ nday }.cueLoc( UEs{ nday }.isHoldTrial);
end

%%
SmoothedSpikingMatrix = AllTypeSmoothedSpikingMatrix(alignInds,:,:);
WholeLFADSRatesMatrix = AllTypeWholeLFADSRatesMatrix(alignInds,:,:);
nTrials = length(clVector);
nNeurons = size(alf{nday}(1).binned_spikes, 1);
nTimesLFADS = size(alf{nday}(1).rates, 2);
%%
AllCueLoc = unique(clVector);
%rebinSize = par.spikeBinMs; % data was rebinned from 1ms to 10ms during LFADS
AvgFiringRate = zeros( size( alf{ nday }( 1 ).binned_spikes ) );
% initialize a avg firing rate matrix to store avg true firing rate (nNeurons x n Rebinned Times)
% for n = 1:nNeuron % loop over all neurons
AvgLFADSRate = zeros( size( alf{ nday }( 1 ).rates ) );
% initialize a avg LFADS rate matrix to store avg LFADS rates (nNeurons x n Rebinned Times)
AvgFiringRate_InRF = zeros(size(AvgFiringRate));
AvgFiringRate_OffRF = zeros(size(AvgFiringRate));

% initialize avg Firing rates matrix for differet cue location
AvgLFADSRate_InRF = zeros(size(AvgLFADSRate));
AvgLFADSRate_OffRF = zeros(size(AvgLFADSRate));

% initialize avg LFADS rates matrix for differet cue location
TrialIndexInRF = false(nNeurons, nTrials);
TrialIndexOffRF = false(nNeurons, nTrials);
% initialize a matrix for each cue location to store the trial index. Ones
% would be index for trials that have the corresponding cue location

for n = 1:nNeurons % loop over all neurons

    AvgFiringRate(n,:) = (sum(SmoothedSpikingMatrix(:,:,n),1))*(1/nTrials)*(1000/binsize_rescaled);
    % Store the avg firing rate to the avgFiringRate matrix
    AvgLFADSRate(n,:) = (sum(WholeLFADSRatesMatrix(:,:,n),1))*(1/nTrials);

    TrialIndexInRF(n,:) = clVector == rfLoc{nday}(n,1);
    nInRF = length(clVector(clVector == rfLoc{nday}(n,1)));
    % find the number of trials that cue loc is in RF
    TrialIndexOffRF(n,:) = clVector == rfLoc{nday}(n,2);
    % Find the index of the trials that cue loc is in Rf or off Rf
    nOffRF = length(clVector(clVector == rfLoc{nday}(n,2)));
    % find the number of trials that cue loc is Off RF
%     OutRF = AllCueLoc(~ismember(AllCueLoc, rfLoc(n,:))); 
%     % Find the cue locations that are not in or off Rf of this neuron
%     TrialIndexOutRF1(n,:) = clVector == OutRF(1);
%     % always find the index of the trials that cue loc number is smaller and
%     % put these index into the OutRF1
%     nOutRF1 = length(clVector(clVector == OutRF(1)));
%     % find the number of trials that cue loc is OutRF 1
%     TrialIndexOutRF2(n,:) = clVector == OutRF(2);
%     % always find the index of the trials that cue loc number is larger and
%     % put these index into the OutRF2
%     nOutRF2 = length(clVector(clVector == OutRF(2)));
%     % find the number of trials that cue loc is OutRF 2

    AvgFiringRate_InRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF)*(1000/binsize_rescaled);
    AvgFiringRate_OffRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF)*(1000/binsize_rescaled);
%     AvgFiringRate_OutRF1(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1)*(1000/par.spikeBinMs);
%     AvgFiringRate_OutRF2(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2)*(1000/par.spikeBinMs);

    AvgLFADSRate_InRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF);
    AvgLFADSRate_OffRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF);
%     AvgLFADSRate_OutRF1(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1);
%     AvgLFADSRate_OutRF2(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2);
end



%% plots of avg neuron firing rate and avg LFADS rate for different cue locations， plus spiking rasters
%savedirOne = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDa%y_CO_AO_TD_HoldRel_JanToApr/' ...
%    'postAnalysis/withGoodNeurons_PBTRun_20180403/PSTH/NoSmoothing/cueOnset/170127/'];
savedirOne = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withGoodNeurons_PBTRun_20180403/PSTH/committeeMeetingVersion';

cd(savedirOne);
clear set

for n = 1:nNeurons
    f1 = figure
    y_max1 = max([AvgFiringRate_InRF(n,:), AvgFiringRate_OffRF(n,:)])
    y_max2 = max([AvgLFADSRate_InRF(n,:), AvgLFADSRate_OffRF(n,:)])
    y_max = max(y_max1, y_max2)
    
    sPlot_avgRealRates = subplot(3,2,1) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
%     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
%     hold on
%     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-800','cueOnset','+800'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgRealRates, 'Avg Real Firing Rate', 'FontSize', 16);
    ylabel('Firing Rate');
    axes(sPlot_avgRealRates)
    ylim([0 y_max])
    
%     xlabel('time (ms)');
%     legend('show')
    
    sPlot_avgLFADSRates = subplot(3,2,2) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgLFADSRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgLFADSRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
%     plot(AvgLFADSRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
%     hold on
%     plot(AvgLFADSRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS])
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','cueOnset','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgLFADSRates, 'Avg LFADS Firing Rate', 'FontSize', 16);
    ylabel('Firing Rate');
    %     set(gca,'YLim',sPlot_avgRealRates.YLim);
    axes(sPlot_avgLFADSRates)
    ylim([0 y_max])
    xlabel('time (ms)');
    Legend = legend('show');
    set(Legend,'Position',[0.7802 0.82 0.1226 0.0425]);
    
    
    sPlot_InRFRaster = subplot(3,2,3)
    imagesc(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-800','cueOnset','+800'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    sPlot_LFADSInRFRaster = subplot(3,2,4)
    imagesc(WholeLFADSRatesMatrix(TrialIndexInRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-800','cueOnset','+800'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_LFADSInRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    sPlot_OffRFRaster = subplot(3,2,5)
    imagesc(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-800','cueOnset','+800'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
    sPlot_LFADSOffRFRaster = subplot(3,2,6)
    imagesc(WholeLFADSRatesMatrix(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-800','cueOnset','+800'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_LFADSOffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
   
    
%     sPlot_OutRF1Raster = subplot(5,2,7)
%     imagesc(SmoothedSpikingMatrix(TrialIndexOutRF1(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
%     % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
% %     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
%     set(gca,'xticklabel',{[]})
%     % Put the appropriate labels on the x axis. * remember to change this
%     % if necessary
%     title(sPlot_OutRF1Raster, 'Cue outside RF 1');
%     ylabel('Trials');
%     % ylim([0 nOutRF1]);
%     
%     
%     sPlot_OutRF2Raster = subplot(5,2,9)
%     imagesc(SmoothedSpikingMatrix(TrialIndexOutRF2(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
%     % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-750','CueOnset','+750'});
%     % Put the appropriate labels on the x axis. * remember to change this
%     % if necessary
%     title(sPlot_OutRF2Raster, 'Cue outside RF 2');
%     ylabel('Trials');
%     % ylim([0 nOutRF2]);
%     xlabel('time (ms)');
    
%     suptitle(['Multi-unit ' char(neuronID(n))]);
%     set(f1, 'Position', [428 242 1215 811]);
%     
%     print(f1,['Neuron ' char(neuronID(n))], '-dpng');
    suptitle(['Multi-unit ' int2str(n)]);
    set(f1, 'Position', [428 142 1215 811]);
    print(f1,['Multi-unit ' int2str(n)], '-dpng');

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
    set(f2, 'Position', [297 203 1315 796]);
    
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