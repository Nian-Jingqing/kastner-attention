function  plottingPSTHs_fourLocations(alignedStruct, nTimesLFADS, alignedLabel, nIndices, savedirOne  )
%PLOTTINGPSTHS Summary of this function goes here
%   Detailed explanation goes here
    realRate_cond1 = alignedStruct.realRate_cond1;
    realRate_cond2 = alignedStruct.realRate_cond2;
    realRate_cond3 = alignedStruct.realRate_cond3;
    realRate_cond4 = alignedStruct.realRate_cond4;
    realRaster = alignedStruct.realRaster;
    lfadsRate_cond1 = alignedStruct.lfadsRate_cond1;
    lfadsRate_cond2 = alignedStruct.lfadsRate_cond2;
    lfadsRate_cond3 = alignedStruct.lfadsRate_cond3;
    lfadsRate_cond4 = alignedStruct.lfadsRate_cond4;
    lfadsRaster = alignedStruct.lfadsRaster;
    TrialIndexInRF = alignedStruct.trialIndexInRF;
    TrialIndexOffRF = alignedStruct.trialIndexOffRF;
    TrialIndexPlusRF = alignedStruct.trialIndexPlusRF;
    TrialIndexMinusRF = alignedStruct.trialIndexMinusRF;

    
nNeurons = size(realRate_cond1, 1);
% for n = 3
for n = 1:nNeurons
    f1 = figure;
    % set s1 and s2 to have same y limit
    y_max1 = max([realRate_cond1(n,:), realRate_cond2(n,:), realRate_cond3(n,:), realRate_cond4(n,:)])
    y_max2 = max([lfadsRate_cond1(n,:), lfadsRate_cond2(n,:), lfadsRate_cond3(n,:), lfadsRate_cond4(n,:)])
    y_max = max(y_max1, y_max2)

    
    s1 = subplot(5,2,1);
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(realRate_cond1(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(realRate_cond2(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on
    plot(realRate_cond3(n,:), 'g', 'DisplayName','Cue at RF plus');
    hold on
    plot(realRate_cond4(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue at RF minus');
    hold on
%     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
%     hold on
%     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s1, 'Avg Real Firing Rate');
    ylabel('Firing Rate');
    axes(s1)
    ylim([0 y_max])
    
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


    s2 = subplot(5,2,2);
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(lfadsRate_cond1(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(lfadsRate_cond2(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on
    plot(lfadsRate_cond3(n,:), 'g', 'DisplayName','Cue at RF plus');
    hold on
    plot(lfadsRate_cond4(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue at RF minus');
    hold on
        %     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
        %     hold on
        %     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
        % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
        %     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
        % Put the appropriate labels on the x axis. * remember to change this
        % if necessary
        %set(gca,'YLim',s1.YLim)
    axes(s2)
    ylim([0 y_max])
    title(s2, 'Avg LFADS Rate');
    ylabel('Firing Rate');

    %     xlabel('time (ms)');
    legend('show')
    set(legend, 'Location', 'northoutside');
                  
    
    s3 = subplot(5,2,3);
    imagesc(realRaster(TrialIndexInRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s3, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    s4 = subplot(5,2,4);
    imagesc(lfadsRaster(TrialIndexInRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s4, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    s5 = subplot(5,2,5);
    imagesc(realRaster(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s5, 'Cue opposite RF');
    ylabel('Trials');
    % xlabel('time (ms)');
    % ylim([0 nInRF]);
    
    
    
    s6 = subplot(5,2,6);
    imagesc(lfadsRaster(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s6, 'Cue opposite RF');
    ylabel('Trials');
    % xlabel('time (ms)');
    % ylim([0 nOffRF]);

    s7 = subplot(5,2,7);
    imagesc(realRaster(TrialIndexPlusRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s7, 'Cue at RF plus');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    s8 = subplot(5,2,8);
    imagesc(lfadsRaster(TrialIndexPlusRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s8, 'Cue at RF plus');
    ylabel('Trials');
    % ylim([0 nInRF]);

    s9 = subplot(5,2,9);
    imagesc(realRaster(TrialIndexMinusRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s9, 'Cue at RF minus');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nInRF]);
    
    
    
    s10 = subplot(5,2,10);
    imagesc(lfadsRaster(TrialIndexMinusRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{'-300',alignedLabel,'+300'});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s10, 'Cue at RF minus');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
   
    suptitle(['Multi-unit ' int2str(nIndices(n))]);
    set(f1, 'Position', [310 4 1215 962]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');

    close;
end




end