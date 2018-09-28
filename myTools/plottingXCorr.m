function plottingXCorr(n, xcorr_struct, condType, color, infield, titleKeyword, timeLag)
    timeLagStr = num2str(timeLag);
    addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
    shuffledFieldName = [infield, '_shuffled'];
    shadedErrorBar([], xcorr_struct(n).(infield), {@mean, @(x) std(x)./sqrt(size(xcorr_struct(n).(infield), 1)) }, 'lineProps', color);
    shadedErrorBar([], xcorr_struct(n).(shuffledFieldName), {@mean, @(x) std(x)./sqrt(size(xcorr_struct(n).(shuffledFieldName), 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});
    set(gca,'XTick',[1 0.5*size(xcorr_struct(n).(infield), 2) size(xcorr_struct(n).(infield), 2)]);
    timeStrLeft = ['-', timeLagStr];
    timeStrRight = ['+', timeLagStr];
    set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
    set(gca,'XLim',[0 size(xcorr_struct(n).(infield), 2)])
    ylabel('cross-correlation');
    title(['Cross-correlation Between ' titleKeyword ' and LFP for ' condType])
end