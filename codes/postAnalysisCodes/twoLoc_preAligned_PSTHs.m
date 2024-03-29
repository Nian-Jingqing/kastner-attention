if alignType == 1
    alignTime = cueOnset;
    alignStr = 'cueOnset';
    leftBound = -300;
    rightBound = +900;
    savedir1_base = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20181023/PSTH/cueOnset/';
elseif alignType == 2
    alignTime = arrayOnset;
    alignStr = 'arrayOnset';
    leftBound = -400;
    rightBound = +800;
    savedir1_base = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20181023/PSTH/arrayOnset/';
else
    alignTime = dim;
    alignStr = 'targetDim';
    leftBound = -600;
    rightBound = +600;
    savedir1_base = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20181023/PSTH/targetDim/';
end

%%
% for nday = 1: numel(datasets)
getGoodNeurons
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
if alignType == 1
    alignInds = ( alignType - 1 ) * numel( UEs{ nday }.cueLoc ) + (1 : numel( UEs{ nday }.cueLoc ));
    clVector = UEs{ nday }.cueLoc;
elseif alignType == 3
    alignInds = ( alignType - 1 ) * numel( UEs{ nday }.cueLoc ) + (1 : numel( find(UEs{ nday }.isHoldTrial )));
    clVector = UEs{ nday }.cueLoc( UEs{ nday }.isHoldTrial);
else
    alignInds = ( alignType - 1 ) * numel( UEs{ nday }.cueLoc ) + find(UEs{nday}.isHoldTrial);
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
savedir1 = [savedir1_base datasets( nday ).shortName];
cd(savedir1);
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
    set(gca,'XTick',[1 alignTime nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{int2str(leftBound), alignStr,int2str(rightBound)});
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
    set(gca,'XTick',[1 alignTime nTimesLFADS])
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{int2str(leftBound), alignStr,int2str(rightBound)}); 
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
    set(gca,'XTick',[1 alignTime nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{int2str(leftBound), alignStr,int2str(rightBound)});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    sPlot_LFADSInRFRaster = subplot(3,2,4)
    imagesc(WholeLFADSRatesMatrix(TrialIndexInRF(n,:),:,n)) 
    set(gca,'XTick',[1 alignTime nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{int2str(leftBound), alignStr,int2str(rightBound)});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_LFADSInRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    sPlot_OffRFRaster = subplot(3,2,5)
    imagesc(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 alignTime nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{int2str(leftBound), alignStr,int2str(rightBound)});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
    sPlot_LFADSOffRFRaster = subplot(3,2,6)
    imagesc(WholeLFADSRatesMatrix(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 alignTime nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{int2str(leftBound), alignStr,int2str(rightBound)});
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
    suptitle(['Multi-unit ' int2str(nIndices(n))]);
    set(f1, 'Position', [428 142 1215 811]);
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');

    close;
end


%%