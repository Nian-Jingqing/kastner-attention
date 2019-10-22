function plottingSpikePanel( bb, fig_title, date )

%% Remove NaN for minSpikeBand
chVsb{32} = 0;
for ich = 1:size(bb.validSpikeBand,2)
    tic;
    tmp = bb.validSpikeBand(9062009 : 18124016, ich);
    %whereNan_vsb = find( isnan( bb.validSpikeBand( : , ich ) ) );
    whereNan_vsb = find( isnan( tmp ) );
    %chVsb{ ich } = bb.validSpikeBand(1:(whereNan_vsb(1) - 1), ich);
    if ~isempty(whereNan_vsb)
        chVsb{ ich } = tmp(1:(whereNan_vsb(1) - 1));
    else
        chVsb{ ich } = tmp;
    end
    toc;
end
%%
vsb_mat = [chVsb{:}];

%% Call CP's function to grab threshold crossings (NEED TO WORK ON THIS)

% should call the function calcThresholdCrossings. The function takes N x T data, but my data is in cell array. Need to check whether all channels have same length.
threshMultOrFixed = -5;
useMultiplier = true;
tStep = 1/40000;
windowLength = 40; % default is 30. Changed to 40 to match a 40000-equvalent sampling rate as default
wf = vsb_mat';

[t,inds, wfstds, threshMultOrFixed] = calcThresholdCrossings(wf, threshMultOrFixed, windowLength, tStep, useMultiplier);

%% make plots (NEED TO WORK ON THIS)
figure
for ich = 1:32
    subplot(4,8,ich)
    for iSpike = 1:numel(inds{ich})
        plot(vsb_mat((inds{ich}(iSpike)-10:inds{ich}(iSpike)+10), ich)')
        hold on
    end
    title(['Channel ' int2str(ich)])
    axis tight
end
suptitle(fig_title)
set(gcf, 'Position', [4 4 1914 1082])

savedir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/';
if ~isdir(savedir)
    mkdir(savedir);
end
cd(savedir)
saveFileName = ['Spiking_panel_', date];
print(gcf, saveFileName, '-dpng');