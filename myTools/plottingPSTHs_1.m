function  plottingPSTHs_1(realRate_cond1,realRate_cond2, realRaster, lfadsRate_cond1,lfadsRate_cond2, lfadsRaster, TrialIndexInRF, TrialIndexOffRF, nTimesLFADS, alignedLabel, nIndices, savedirOne, r, left, right)
%PLOTTINGPSTHS Summary of this function goes here
%   Detailed explanation goes here

nNeurons = size(realRate_cond1, 1);
% for n = 3
for n = 1:nNeurons
%for n = 23
    f1 = figure;
    y_max1 = max([realRate_cond1(n,:), realRate_cond2(n,:)])
    y_max2 = max([lfadsRate_cond1(n,:), lfadsRate_cond2(n,:)])
    y_max = max(y_max1, y_max2)
    %y_max = 88;
    s1 = subplot(3,2,1);
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(realRate_cond1(n,:), 'r', 'LineWidth', 2, 'DisplayName','Cue in RF');
    hold on
    plot(realRate_cond2(n,:), 'b', 'LineWidth', 2, 'DisplayName','Cue opposite RF');
    hold on 
%     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
%     hold on
%     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 r*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{left,alignedLabel,right});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s1, 'Avg Real Firing Rate', 'FontSize', 16);
    ylabel('Firing Rate');
    axes(s1)
    ylim([0 y_max])
    set(gca, 'fontsize', 12);    
    % xlabel('time (ms)');
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


    s2 = subplot(3,2,2);
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(lfadsRate_cond1(n,:), 'r', 'LineWidth', 2, 'DisplayName','Cue in RF');
    hold on
    plot(lfadsRate_cond2(n,:), 'b', 'LineWidth', 2, 'DisplayName','Cue opposite RF');
    hold on 
        %     plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
        %     hold on
        %     plot(AvgFiringRate_OutRF2(n,:), 'Color', [1, 0.8, 0], 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 r*nTimesLFADS nTimesLFADS]);
        % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
        %     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{left,alignedLabel,right});
        % Put the appropriate labels on the x axis. * remember to change this
        % if necessary
    set(gca,'YLim',s1.YLim)
    title(s2, 'Avg LFADS Rate', 'FontSize', 16);
    ylabel('Firing Rate');
    axes(s2)
    ylim([0 y_max])
    %     xlabel('time (ms)');
    set(gca, 'fontsize', 12);    
    %legend('show')
    %set(legend, 'Location', 'best'); 
    %set(legend, 'Location', 'southeast');
    
    


    
    s3 = subplot(3,2,3);
    %c = hot(20);
    imagesc(realRaster(TrialIndexInRF(n,:),:,n))
    c = flipud(gray);
    c_1 = [c(1:3,:); c(38:64,:)];
    colormap(gca, c_1)
    %colormap(flipud(gray))
    %colormap(c);
    %cmap = colormap;
    %cmap_1 = [cmap(1:3,:); cmap(50:64,:)];
    %colormap(cmap_1);
    %keyboard
    set(gca,'XTick',[1 r*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{left,alignedLabel,right});
    set(gca, 'fontsize', 12);    
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s3, 'Cue in RF');
    ylabel('Trials');
    % xlabel('time (ms)');
    % ylim([0 nInRF]);
    
    
    s4 = subplot(3,2,4);
    imagesc(lfadsRaster(TrialIndexInRF(n,:),:,n))
    %c_2 = [c(1:25,:); c(35:64,:)];
    %colormap(gca, c_2)
    colormap(gca, flipud(gray))
    set(gca,'XTick',[1 r*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{left,alignedLabel,right});
    set(gca, 'fontsize', 12);    
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s4, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    
    s5 = subplot(3,2,5);
    imagesc(realRaster(TrialIndexOffRF(n,:),:,n))
    colormap(gca, c_1)
    set(gca,'XTick',[1 r*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{left,alignedLabel,right});
    set(gca, 'fontsize', 12);    
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s5, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nInRF]);
    
    
    
    s6 = subplot(3,2,6);
    imagesc(lfadsRaster(TrialIndexOffRF(n,:),:,n))
    colormap(gca, flipud(gray))
    set(gca,'XTick',[1 r*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'XTickLabels',{left,alignedLabel,right});
    set(gca, 'fontsize', 12);
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(s6, 'Cue opposite RF');
    ylabel('Trials');
    xlabel('time (ms)');
    % ylim([0 nOffRF]);
    
    
   
    suptitle(['Multi-unit ' int2str(nIndices(n))]);
    set(f1, 'Position', [428 200 1215 811]);
    cd(savedirOne)
    print(f1,['Multi-unit ' int2str(nIndices(n))], '-dpng');
    printpdf(f1,int2str(nIndices(n)) )
    close;
end




end

