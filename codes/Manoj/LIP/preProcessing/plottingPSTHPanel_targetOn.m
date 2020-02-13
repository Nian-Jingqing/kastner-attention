function plottingPSTHPanel_targetOn(bb, date, multiple)
UE_baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/';
UE_file = ['UE_' date '.mat'];
load_UE_path = fullfile(UE_baseDir, UE_file);
clear UE
%load_UE_path = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/UE_03062019.mat';
load(load_UE_path);

startInds = UE.fixOn;
stopInds = UE.trialEnd + 400;
trialstruct = struct;
for itrial = 1:numel(startInds)
    trialstruct(itrial).isErrorTrial = UE.isErrorTrial(itrial);
    trialstruct(itrial).isEarlyError = UE.isEarlyError(itrial);
    trialstruct(itrial).cueType = UE.cueType(itrial);
    trialstruct(itrial).barType = UE.cueType(itrial);
    trialstruct(itrial).fixType = UE.fixType(itrial);
    trialstruct(itrial).isValidTarget = UE.isValidTarget(itrial);
    trialstruct(itrial).isSameObjTarget = UE.isSameObjTarget(itrial);
    trialstruct(itrial).isDiffObjTarget = UE.isDiffObjTarget(itrial);    
    trialstruct(itrial).startTime = startInds(itrial);
    trialstruct(itrial).endTime = stopInds(itrial);
    trialstruct(itrial).startInd = startInds(itrial);
    trialstruct(itrial).endInd = stopInds(itrial);
    trialstruct(itrial).condition = UE.cueType(itrial); % need to change!!        
end


%% get threshold to use for each channel
thresh = zeros(1, size(bb.minSpikeBand, 2)); % for each channel, there is one threshold
seg_len = 9062008;
num_segs = 5;
for iseg = 1:num_segs
    startInd = ((iseg - 1) * seg_len + 1);
    endInd = startInd + seg_len - 1;
    tmp = bb.validSpikeBand(startInd : endInd, :);
    wfstds = std(tmp, 0, 1);
    thresh = thresh + wfstds;
    clear tmp
    iseg
end

thresh = (-1) * multiple * (thresh / num_segs);

%% get var and remove NaN, also remove NaN for minSpikeBand
for ich = 1:size(bb.minSpikeBand,2)
    chVar{ ich } = bb.meanSquared( bb.meanSquaredChannel == ich );
    whereNan = find( isnan( chVar{ ich } ) );
    chVar{ ich } = chVar{ ich }( 1 : (whereNan(1) - 1) );
    %    whereNan_msb = find( isnan( bb.minSpikeBand_rmCA( : , ich ) ) );
    whereNan_msb = find( isnan( bb.minSpikeBand( : , ich ) ) );
    %chMsb{ ich } = bb.minSpikeBand_rmCA(1:(whereNan_msb(1) - 1), ich);
    chMsb{ ich } = bb.minSpikeBand(1:(whereNan_msb(1) - 1), ich);
end


%% get spikes
spikes = sparse(size(chMsb{1}, 1), numel(chMsb));
for ich = 1:size(bb.minSpikeBand,2)
    chThres = thresh(ich);
    leftValue = chMsb{ich} - chThres;
    spikes(leftValue <= 0, ich) = 1;
end
% remove spikes if show on more than 25 TCs (~80%)
allSpikesMS = sum(spikes,2);
spikes((allSpikesMS > 25), :) = 0;
% % get rid of MU 1 - 4
%spikes = spikes(:, 5:end);
stream.spikes = sparse(spikes);

%%
dtMS = 1;
C = Continuous.Continuous(stream, dtMS);
sigma_neural = 10;
C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );
r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );
tc_r = R.Rstruct(r.r);
%
sigma = 100;
binsize_rescaled = 10;
%
tc_r.smoothFieldInR('spikes', 'spike_smoothed', sigma, 1);
tc_rebinned = tc_r.binData({'spike_smoothed'}, [binsize_rescaled]);

barStart = UE.barOn - UE.fixOn;
cueStart = UE.cueOn - UE.fixOn;
targetStart = UE.targetOn - UE.fixOn;
tc = [];
window_preCue = round([-300  0] / binsize_rescaled);
timePoints_preCue = window_preCue(1):window_preCue(2);
preCueTensor = zeros(numel(tc_rebinned), size(tc_rebinned(1).spike_smoothed, 1), numel(timePoints_preCue));
keyboard
for ntr = 1:numel(tc_rebinned)
    tc(ntr).spikes = tc_rebinned(ntr).spike_smoothed;
    tc( ntr ).barOnset = round( barStart( ntr ) /binsize_rescaled);
    tc( ntr ).cueOnset = round( cueStart( ntr ) / binsize_rescaled);
    tc( ntr ).targetStart = round( targetStart( ntr ) / binsize_rescaled);
    preCueTensor(ntr, :,:) = tc(ntr).spikes(:, tc(ntr).cueOnset + timePoints_preCue);
end

%$ get PSTHs and rasters
events = {'barOn', 'cueOn', 'targetOn'};
conditions = struct;
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
window = round([-300  500] / binsize_rescaled);
timePoints = window(1):window(2);

%
for nBar = 1:numel(conditions.barOn)
    trialIndicesPerCond.barOn.(conditions.barOn{nBar}) = UE.barType == nBar;
    % compute tc psth
    [psth, raster_tensor, stderr] = preparePSTH_Manoj(tc, 'spikes', trialIndicesPerCond.barOn.(conditions.barOn{nBar}), [tc.barOnset], timePoints, binsize_rescaled);
    PSTH.barOn.(conditions.barOn{nBar}).tc.raster_tensor = raster_tensor;
    PSTH.barOn.(conditions.barOn{nBar}).tc.psth = psth;
    PSTH.barOn.(conditions.barOn{nBar}).tc.stderr = stderr;
    
    for nCue = 1 : numel(conditions.cueOn)
        cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
        trialIndicesPerCond.cueOn.(cue_cond_field) = (UE.barType == nBar) & (UE.cueType == nCue);

        % compute lfads psth
        [psth, raster_tensor, stderr] = preparePSTH_Manoj(tc, 'spikes', trialIndicesPerCond.cueOn.(cue_cond_field), [tc.cueOnset], timePoints, binsize_rescaled);
        PSTH.cueOn.(cue_cond_field).tc.raster_tensor = raster_tensor;
        PSTH.cueOn.(cue_cond_field).tc.psth = psth;
        PSTH.cueOn.(cue_cond_field).tc.stderr = stderr;
    end
end

%% plotting PSTHs
savedir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/newSessions/';
if ~isdir(savedir)
    mkdir(savedir);
end

pre_post_times = [300, 500] / binsize_rescaled;
t = pre_post_times;
t = [-t(1):-1 0:t(2)] * binsize_rescaled; %t(end)=[];

whichEvent = 'cueOn';
cond_1 = 'Exo_TL_Vert';
cond_2 = 'Exo_TR_Vert';
cond_3 = 'Exo_BL_Vert';
cond_4 = 'Exo_BR_Vert';

%keyboard

%%
figure
for ich = 1:32
    subplot(4,8,ich)
    p_1 = PSTH.(whichEvent).(cond_1).tc.psth(ich, :);
    p_1_upper = p_1 + PSTH.(whichEvent).(cond_1).tc.stderr(ich, :);
    p_1_lower = p_1 - PSTH.(whichEvent).(cond_1).tc.stderr(ich, :);
    p_2 = PSTH.(whichEvent).(cond_2).tc.psth(ich, :);
    p_2_upper = p_2 + PSTH.(whichEvent).(cond_2).tc.stderr(ich, :);
    p_2_lower = p_2 - PSTH.(whichEvent).(cond_2).tc.stderr(ich, :);
    p_3 = PSTH.(whichEvent).(cond_3).tc.psth(ich, :);
    p_3_upper = p_3 + PSTH.(whichEvent).(cond_3).tc.stderr(ich, :);
    p_3_lower = p_3 - PSTH.(whichEvent).(cond_3).tc.stderr(ich, :);
    p_4 = PSTH.(whichEvent).(cond_4).tc.psth(ich, :);
    p_4_upper = p_4 + PSTH.(whichEvent).(cond_4).tc.stderr(ich, :);
    p_4_lower = p_4 - PSTH.(whichEvent).(cond_4).tc.stderr(ich, :);    
    
    h(1) = plot(t, p_1, '-b', 'LineWidth', 2, 'DisplayName', 'Exo-TL-Vert')
    hold on
    plot(t, p_1_upper, '-b', 'LineWidth', 1)
    hold on
    plot(t, p_1_lower, '-b', 'LineWidth', 1)
    hold on
    h(2) = plot(t, p_2, '-r', 'LineWidth', 2, 'DisplayName', 'Exo-TR-Vert')
    hold on
    plot(t, p_2_upper, '-r', 'LineWidth', 1)
    hold on
    plot(t, p_2_lower, '-r', 'LineWidth', 1)
    hold on
    h(3) = plot(t, p_3, '-g', 'LineWidth', 2, 'DisplayName', 'Exo-BL-Vert')
    hold on
    plot(t, p_3_upper, '-g', 'LineWidth', 1)
    hold on
    plot(t, p_3_lower, '-g', 'LineWidth', 1)
    hold on
    h(4) = plot(t, p_4, '-y', 'LineWidth', 2, 'DisplayName', 'Exo-BR-Vert')
    hold on
    plot(t, p_4_upper, '-y', 'LineWidth', 1)
    hold on
    plot(t, p_4_lower, '-y', 'LineWidth', 1)
    hold on    
    
    yl = ylim;
    plot([0 0], yl, '--k')
    if ich == 1
        %legend(cond_1, cond_2, whichEvent)
        %legend('show')
        legend('Location', 'best')
        legend(h(1:4));
    end
    title(['Channel ' int2str(ich)])
    xlabel('(ms)')
    axis tight
end
suptitle('Exo-Vertical Bars')
set(gcf, 'Position', [4 4 1914 1082]);
%tightfig(gcf)
cd(savedir);
a = sprintf('%.2f', multiple);
saveFileName = [date, '_', a, 'std_PSTH_panel_', 'Exo_Vert'];
print(gcf, saveFileName, '-dpng');

%%
clear h
cond_1 = 'Exo_TL_Hori';
cond_2 = 'Exo_TR_Hori';
cond_3 = 'Exo_BL_Hori';
cond_4 = 'Exo_BR_Hori';
figure
for ich = 1:32
    subplot(4,8,ich)
    p_1 = PSTH.(whichEvent).(cond_1).tc.psth(ich, :);
    p_1_upper = p_1 + PSTH.(whichEvent).(cond_1).tc.stderr(ich, :);
    p_1_lower = p_1 - PSTH.(whichEvent).(cond_1).tc.stderr(ich, :);
    p_2 = PSTH.(whichEvent).(cond_2).tc.psth(ich, :);
    p_2_upper = p_2 + PSTH.(whichEvent).(cond_2).tc.stderr(ich, :);
    p_2_lower = p_2 - PSTH.(whichEvent).(cond_2).tc.stderr(ich, :);
    p_3 = PSTH.(whichEvent).(cond_3).tc.psth(ich, :);
    p_3_upper = p_3 + PSTH.(whichEvent).(cond_3).tc.stderr(ich, :);
    p_3_lower = p_3 - PSTH.(whichEvent).(cond_3).tc.stderr(ich, :);
    p_4 = PSTH.(whichEvent).(cond_4).tc.psth(ich, :);
    p_4_upper = p_4 + PSTH.(whichEvent).(cond_4).tc.stderr(ich, :);
    p_4_lower = p_4 - PSTH.(whichEvent).(cond_4).tc.stderr(ich, :);    
    
    h(1) = plot(t, p_1, '-b', 'LineWidth', 2, 'DisplayName', 'Exo-TL-Hori')
    hold on
    plot(t, p_1_upper, '-b', 'LineWidth', 1)
    hold on
    plot(t, p_1_lower, '-b', 'LineWidth', 1)
    hold on
    h(2) = plot(t, p_2, '-r', 'LineWidth', 2, 'DisplayName', 'Exo-TR-Hori')
    hold on
    plot(t, p_2_upper, '-r', 'LineWidth', 1)
    hold on
    plot(t, p_2_lower, '-r', 'LineWidth', 1)
    hold on
    h(3) = plot(t, p_3, '-g', 'LineWidth', 2, 'DisplayName', 'Exo-BL-Hori')
    hold on
    plot(t, p_3_upper, '-g', 'LineWidth', 1)
    hold on
    plot(t, p_3_lower, '-g', 'LineWidth', 1)
    hold on
    h(4) = plot(t, p_4, '-y', 'LineWidth', 2, 'DisplayName', 'Exo-BR-Hori')
    hold on
    plot(t, p_4_upper, '-y', 'LineWidth', 1)
    hold on
    plot(t, p_4_lower, '-y', 'LineWidth', 1)
    hold on    
    
    yl = ylim;
    plot([0 0], yl, '--k')
    if ich == 1
        %legend(cond_1, cond_2, whichEvent)
        %legend('show')
        legend(h(1:4));
        legend('Location', 'best')
    end
    title(['Channel ' int2str(ich)])
    xlabel('(ms)')
    axis tight
end
set(gcf, 'Position', [4 4 1914 1082]);
tightfig(gcf)
suptitle('Exo-Horizontal Bars')
cd(savedir);
a = sprintf('%.2f', multiple);
saveFileName = [date, '_', a, 'std_PSTH_panel_', 'Exo_Hori'];
%saveFileName = ['PSTH_panel_', 'Exo_Hori_', date];
print(gcf, saveFileName, '-dpng');
clear all