%% add dataset path
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')

%%
loadpath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
    'multi-unit/continuousOverlapChop/multiDay_JanToMar/updated_Rstruct/170127_cueOnArrayOnTargetDim_HoldRel.mat'];

olapChopped = load(loadpath);
olapChopped = olapChopped.combinedData;
%%
trueSpikes = [olapChopped.r.r.spikes];
nSpikes = size(trueSpikes, 2);
trial_time_ms = 500;
trial_olap_ms = 100;
out = olapChopped.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );

%%
tot_spikes = out.counts;
% nPieces = size(tot_spikes, 1);
%spikes = zeros(size(trueSpikes, 1), nSpikes);
spikes(length(olapChopped.r.r)).spikes = 1;
%% assembly
nChops = 0;
for r_trial = 1: length(olapChopped.r.r)
    spikes(r_trial).spikes = zeros(size(olapChopped.r.r(r_trial).spikes));
    nPieces = ceil( (size( spikes(r_trial).spikes,2 ) - trial_time_ms ) / ( trial_time_ms - trial_olap_ms ) );
    for i = 1:nPieces
        spikesThisPiece = squeeze( tot_spikes( i+nChops, :, : ) );
        trialStart = 1 + (i-1)*trial_time_ms - (i-1)*trial_olap_ms;
        trialEnd = trialStart + trial_time_ms - 1;
%         if i~= size( tot_spikes, 1 )
%             trialStart = 1 + (i-1)*trial_time_ms - (i-1)*trial_olap_ms;
%             trialEnd = trialStart + trial_time_ms - 1;
%         else
%             trialStart = nSpikes - trial_time_ms + 1;
%             trialEnd = nSpikes;
%         end
        spikes(r_trial).spikes( :, trialStart:trialEnd ) = spikesThisPiece;
    end
    nChops = nChops + nPieces;
end

% %% trialize the continuous spikes
% r(length(olapChopped.r.r)).spikes = 1;
% past_nSpikes = 0;
% for i = 1: length(r)
%     nSpikesThisTrial = size(olapChopped.r.r(i).spikes,2);
%     start_ind = past_nSpikes +1;
%     end_ind = past_nSpikes + nSpikesThisTrial;
%     r(i).spikes = spikes(:, start_ind:end_ind);
%     past_nSpikes = end_ind;
% end
%%
assembled = R.Rstruct(spikes);

%% get timing info
loadpath = ['/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/' ...
    'M20170127/Gratings/M20170127_PUL_1M-g3-g4-g5-evokedSpiking-v10.mat'];
data = load(loadpath);
UE = data.UE;
clear data

%% get cueOn, arrayOn, targetDim timing
rawSampleRate = 1000;
sessStartTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1)*rawSampleRate;
% session start time - I set it to be the time when the monkey starts
% fixation in the first trial
sessEndTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
extraEndMs = 500;

startInds = round(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue * rawSampleRate - sessStartTime);
startInds(1) = 1;
stopInds = round(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice * rawSampleRate - sessStartTime);

cueInds = round(UE.cueOnset * rawSampleRate - sessStartTime);
arrayInds = round(UE.arrayOnset * rawSampleRate - sessStartTime);
dimInds = round(UE.targetDim * rawSampleRate - sessStartTime);

cueStart = cueInds - startInds;
arrayStart = arrayInds - startInds;
dimStart = dimInds - startInds(UE.isHoldTrial);

stopTimes = stopInds(UE.isHoldTrial) - startInds(UE.isHoldTrial);
%% set up aligned trials

half_ms = 300;
assembled_hold = R.Rstruct(assembled.r(UE.isHoldTrial));
assembled_hold = R.Rstruct(assembled_hold.r((stopTimes - dimStart) > 302));
olapChopped_hold = olapChopped.r.r(UE.isHoldTrial);
olapChopped_hold = olapChopped_hold((stopTimes - dimStart) > 302);
dimStart = dimStart((stopTimes - dimStart) > 302);

% for cueAlign and arrayAlign
for i = 1:length(assembled.r)
    assembled.r(i).cueAlign = assembled.r(i).spikes(:, (cueStart(i)-(half_ms-1)):(cueStart(i)+half_ms) );
    assembled.r(i).arrayAlign = assembled.r(i).spikes(:, (arrayStart(i)-(half_ms-1)):(arrayStart(i)+half_ms) );
end

for i = 1:length(assembled_hold.r)
    assembled_hold.r(i).dimAlign = assembled_hold.r(i).spikes(:, (dimStart(i)-(half_ms-1)):(dimStart(i)+half_ms) );
end

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(assembled.r); % get trial number
nTimesRaw = size(assembled.r(1).cueAlign, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(assembled.r(1).cueAlign, 1); % get neuron nubmer



%% compute avg neuron firing rates for different cue locations

sigma = 10;
binsize = 1;
[AvgFiringRate_InRF_cue,AvgFiringRate_OffRF_cue, SmoothedSpikingMatrix_cue, TrialIndexInRF_cue, TrialIndexOffRF_cue, nTimesLFADS] = computeFRForPSTH(assembled, 'cueAlign', sigma, binsize, olapChopped.r.r);
[AvgFiringRate_InRF_array,AvgFiringRate_OffRF_array, SmoothedSpikingMatrix_array, TrialIndexInRF_array, TrialIndexOffRF_array, nTimesLFADS] = computeFRForPSTH(assembled, 'arrayAlign', sigma, binsize, olapChopped.r.r);
[AvgFiringRate_InRF_dim,AvgFiringRate_OffRF_dim, SmoothedSpikingMatrix_dim, TrialIndexInRF_dim, TrialIndexOffRF_dim, nTimesLFADS] = computeFRForPSTH(assembled_hold, 'dimAlign', sigma, binsize, olapChopped_hold);


%% select good neurons (v12) - new good neurons after May 21, 18: why new good neurons???
% still need to update the 170320 - 170407 *****

nIndices = [1 4 6 7 8 11 12 13 14 16 17 19 20 21 23 24 25 26 27 29 30 31 32];
% % 170127
% nIndices = [1 2 4 5 6 7 8 9 12 14 15 16 17 19 20 23 24 25 26 29 30 31 32];
% % 170130
% nIndices = [1 5 7 8 10 12 13 27 30 32];
% % 170201
% nIndices = [3 4 5 7 11 12 15 17 18 23 28 30];
% % 170211
% nIndices = [9 15 18 20 24 28 30 31 32 34 37 38 39 42 43 44 45 46 47 48 52 54 55 58 60 62 63 64];
% % 170308
% nIndices = [2 3 5 8 9 11 13 17 20 22 25 26 32 33 36 38 42 43 50 54 55 56 57 58];
% % 170311
% nIndices = [4 5 9 12 16 20 24 25 29 30 31 34 35 37 38 40 41 44 49 50 51 52 53 54 56 57 59 61];
% % 170320
% nIndices = [1 4 7 8 11 12 13 14 15 16 18 21 25 27 28 29 30 33 34 36 39 42 43 44 45 47 48 49 50 51 52 57 58 60 61 62 63 64];
% % 170324
% nIndices = [1 2 3 4 5 10 11 13 15 16 18 19 21 22 24 25 27 31 32 35 39 40 41 42 43 44 45 47 50 52 53 54 55 56 57 58 59 60 61 62 63 64];
% % 170327
% nIndices = [1 2 3 5 6 7 8 12 14 16 17 18 20 22 25 26 27 28 29 31 32 34 38 39 40 41 47 48 49 50 51 52 53 54 55 56 60];
% % 170329 first half
% nIndices = [1 3 4 5 6 7 8 9 10 11 13 18 19 20 21 22 24 25 27 30 31 32 35 36 38 39 42 44 45 46 47 48 49 51 52 53 54 55 56 57 59 60 61 62 63 64];
% % 170329 second half
% nIndices = [7 9 10 11 13 14 17 24 25 29 30 31 32 34 36 37 41 46 47 50 51 52 55 57 58 59 62 64];
% % 170331 first half
% nIndices = [2 4 8 9 10 13 15 16 18 21 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 49 50 51 52 53 55 56 57 58 60 61 63 64];
% % 170331 second half
% nIndices = [16 19 20 23 24 32];
% % 170407
%% plots of avg neuron firing rate and avg LFADS rate for different cue locations， plus spiking rasters
%savedirOne = ['/snel/share/share/derived/kastner/nonLFADS_analysis/pulvinar/' ...
%     'Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/assembledPSTH/170127/'];


%cd(savedirOne);
clear set
for n = 3
% for n = 1:nNeurons
    f1 = figure
    
    sPlot_avgRealRates = subplot(3,3,1) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF_cue(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF_cue(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
%     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
%     hold on
%     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','cueOnset','+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgRealRates, 'Avg Real Firing Rate');
    ylabel('Firing Rate');
    
%     xlabel('time (ms)');
     %legend('show')
    
%     sPlot_avgLFADSRates = subplot(3,1,2) 
%     % make the first subplot for plotting avg Reak rates for different cue locations
%     plot(AvgLFADSRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
%     hold on
%     plot(AvgLFADSRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
%     hold on 
% %     plot(AvgLFADSRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
% %     hold on
% %     plot(AvgLFADSRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS])
%     % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-300','cueOnset','+300'}); 
%     % Put the appropriate labels on the x axis. * remember to change this
%     % if necessary
%     title(sPlot_avgLFADSRates, 'Avg LFADS Firing Rate');
%     ylabel('Firing Rate');
% %     set(gca,'YLim',sPlot_avgRealRates.YLim);
%     xlabel('time (ms)');
%     Legend = legend('show');
%     set(Legend,'Position',[0.7802 0.82 0.1226 0.0425]);


    sPlot_avgRealRates_array = subplot(3,3,2) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF_array(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF_array(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
        %     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
        %     hold on
        %     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
        % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
        %     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','arrayOnset','+300'});
        % Put the appropriate labels on the x axis. * remember to change this
        % if necessary
    title(sPlot_avgRealRates_array, 'Avg Real Firing Rate');
    ylabel('Firing Rate');

    %     xlabel('time (ms)');
    %legend('show')
    
    
    sPlot_avgRealRates_dim = subplot(3,3,3) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF_dim(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF_dim(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
        %     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
        %     hold on
        %     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
        % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
        %     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','targetDim','+300'});
        % Put the appropriate labels on the x axis. * remember to change this
        % if necessary
    title(sPlot_avgRealRates_dim, 'Avg Real Firing Rate');
    ylabel('Firing Rate');

    %     xlabel('time (ms)');
    legend('show')
    
    


    
    sPlot_InRFRaster = subplot(3,3,4)
    imagesc(SmoothedSpikingMatrix_cue(TrialIndexInRF_cue(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','cueOnset','+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    sPlot_InRFRaster_array = subplot(3,3,5)
    imagesc(SmoothedSpikingMatrix_array(TrialIndexInRF_array(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','arrayOnset','+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster_array, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    sPlot_InRFRaster_dim = subplot(3,3,6)
    imagesc(SmoothedSpikingMatrix_dim(TrialIndexInRF_dim(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','targetDim','+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster_dim, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    
%     sPlot_LFADSInRFRaster = subplot(3,2,4)
%     imagesc(WholeLFADSRatesMatrix(TrialIndexInRF(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
%     % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
% %     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
%     set(gca,'XTickLabels',{'-800','cueOnset','+800'});
%     % Put the appropriate labels on the x axis. * remember to change this
%     % if necessary
%     title(sPlot_LFADSInRFRaster, 'Cue in RF');
%     ylabel('Trials');
%     % ylim([0 nInRF]);
    
    
    sPlot_OffRFRaster = subplot(3,3,7)
    imagesc(SmoothedSpikingMatrix_cue(TrialIndexOffRF_cue(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','cueOnset','+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
    sPlot_OffRFRaster_array = subplot(3,3,8)
    imagesc(SmoothedSpikingMatrix_array(TrialIndexOffRF_array(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','arrayOnset','+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster_array, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
    sPlot_OffRFRaster_dim = subplot(3,3,9)
    imagesc(SmoothedSpikingMatrix_dim(TrialIndexOffRF_dim(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300','targetDim','+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster_dim, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
    
%     sPlot_LFADSOffRFRaster = subplot(3,2,6)
%     imagesc(WholeLFADSRatesMatrix(TrialIndexOffRF(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
%     % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
% %     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
%     set(gca,'XTickLabels',{'-800','cueOnset','+800'});
%     % Put the appropriate labels on the x axis. * remember to change this
%     % if necessary
%     title(sPlot_LFADSOffRFRaster, 'Cue opposite RF');
%     ylabel('Trials');
%     xlabel('time (ms)');
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
    set(f1, 'Position', [428 200 1215 811]);
    %print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');

%     close;
end





    