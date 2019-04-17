clear spikes
%% set up day information
%datasets(1).shortName = '0208';
%datasets(1).longName = '02082019';
%datasets(2).shortName = '0218';
%datasets(2).longName = '02182019';
datasets(1).shortName = '0218';
datasets(1).midName = '021819';
datasets(1).longName = '02182019';
%datasets(4).shortName = '0227';
%datasets(4).longName = '02272019';
%datasets(5).shortName = '0308';
%datasets(5).longName = '03082019';
%datasets(6).shortName = '0310';
%datasets(6).longName = '03102019';
%datasets(7).shortName = '0311';
%datasets(7).longName = '03112019';


%% load spike sorted data and UE
loadpath = ['/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/'];
disp( sprintf( 'loading chopped day %g / %g', 1, numel( datasets ) ) );
fname = sprintf( '%s%s_v1.mat', loadpath, datasets( 1 ).midName );
tmp = load( fname );
olapChopped_SS = tmp.combinedData;
UE = tmp.combinedData.UE;

%% trialize the threshold crossing (TC) data
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

%% put spike train into a Continuous class
dtMS = 1;
C = Continuous.Continuous(stream, dtMS);
sigma_neural = 10;
C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );

%% turn into a trialized (R) struct
r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

%% combine and save
olapChopped_TC.r = r;
for i = 1:numel(olapChopped_TC.r.r)
    olapChopped_TC.r.r(i).externalInputs = olapChopped_SS.r.r(i).externalInputs;
end

%% process the spikes
smoothOrNot = 1;
sigma = 20;
overWriteRealOrNot = 1;
binsize_rescaled = 1;
trial_time_ms = 500;
trial_olap_ms = 100;
ss = [];
barStart = UE.barOn - UE.fixOn;
cueStart = UE.cueOn - UE.fixOn;
targetStart = UE.targetOn - UE.fixOn;
if overWriteRealOrNot
    % get all real data
    trueSpikes = [olapChopped_SS.r.r.spikes];
    nSpikes = size(trueSpikes, 2);
    out = olapChopped_SS.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );
    tot_spikes = out.counts;
    % nPieces = size(tot_spikes, 1);
    %spikes = zeros(size(trueSpikes, 1), nSpikes);
    spikes(length(olapChopped_SS.r.r)).spikes = 1;
    %% assembly
    nChops = 0;
    for r_trial = 1: length(olapChopped_SS.r.r)
        spikes(r_trial).spikes = zeros(size(olapChopped_SS.r.r(r_trial).spikes));
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
    assembled_real2 = R.Rstruct(spikes);
    if smoothOrNot
        assembled_real2.smoothFieldInR('spikes', 'spike_smoothed', sigma, 1);
        rbinned_real = assembled_real2.binData({'spike_smoothed'}, [binsize_rescaled]);
    else
        rbinned_real = assembled_real2.binData({'spikes'}, [binsize_rescaled]);
    end

    for ntr = 1:numel( olapChopped_SS.r.r )
        if smoothOrNot
            ss( ntr ).spikes = rbinned_real( ntr ).spike_smoothed;
        else            
            ss( ntr ).spikes = rbinned_real( ntr ).spikes;
        end
        ss( ntr ).barOnset = round( barStart( ntr ) /binsize_rescaled);
        ss( ntr ).cueOnset = round( cueStart( ntr ) / binsize_rescaled);
        ss( ntr ).targetStart = round( targetStart( ntr ) / binsize_rescaled);
    end
end

%% do the same for the TC data
tc = [];
if overWriteRealOrNot
    % get all real data
    trueSpikes = [olapChopped_TC.r.r.spikes];
    nSpikes = size(trueSpikes, 2);
    out = olapChopped_TC.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );
    tot_spikes = out.counts;
    % nPieces = size(tot_spikes, 1);
    %spikes = zeros(size(trueSpikes, 1), nSpikes);
    spikes(length(olapChopped_TC.r.r)).spikes = 1;
    %% assembly
    nChops = 0;
    for r_trial = 1: length(olapChopped_TC.r.r)
        spikes(r_trial).spikes = zeros(size(olapChopped_TC.r.r(r_trial).spikes));
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
    assembled_real3 = R.Rstruct(spikes);
    if smoothOrNot
        assembled_real3.smoothFieldInR('spikes', 'spike_smoothed', sigma, 1);
        rbinned_real = assembled_real3.binData({'spike_smoothed'}, [binsize_rescaled]);
    else
        rbinned_real = assembled_real3.binData({'spikes'}, [binsize_rescaled]);
    end

    for ntr = 1:numel( olapChopped_TC.r.r )
        if smoothOrNot
            tc( ntr ).spikes = rbinned_real( ntr ).spike_smoothed;
        else            
            tc( ntr ).spikes = rbinned_real( ntr ).spikes;
        end
    end
end

%%
% basic PSTH save root
saveRoot = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/thresCrossings/verifyTcPSTH/0218/';
if ~isdir(saveRoot)
    mkdir(saveRoot);
end

binsize_rescaled = 1;
events = {'barOn', 'cueOn', 'targetOn'};
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
window = round([-400  400] / binsize_rescaled);
timePoints = window(1):window(2);

%%
saveRootDay = fullfile(saveRoot, datasets(1).shortName);
for nBar = 1:numel(conditions.barOn)
    trialIndicesPerCond.barOn.(conditions.barOn{nBar}) = UE.barType == nBar;
    % compute real psth
    [psth, raster_tensor, stderr] = preparePSTH_Manoj(ss, 'spikes', trialIndicesPerCond.barOn.(conditions.barOn{nBar}), [ss.barOnset], timePoints, binsize_rescaled);
    PSTH.barOn.(conditions.barOn{nBar}).ss.raster_tensor = raster_tensor;
    PSTH.barOn.(conditions.barOn{nBar}).ss.psth = psth;
    PSTH.barOn.(conditions.barOn{nBar}).ss.stderr = stderr;

    % compute lfads psth
    [psth, raster_tensor, stderr] = preparePSTH_Manoj(tc, 'spikes', trialIndicesPerCond.barOn.(conditions.barOn{nBar}), [ss.barOnset], timePoints, binsize_rescaled);
    PSTH.barOn.(conditions.barOn{nBar}).tc.raster_tensor = raster_tensor;
    PSTH.barOn.(conditions.barOn{nBar}).tc.psth = psth;
    PSTH.barOn.(conditions.barOn{nBar}).tc.stderr = stderr;
    
    for nCue = 1 : numel(conditions.cueOn)
        cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
        trialIndicesPerCond.cueOn.(cue_cond_field) = (UE.barType == nBar) & (UE.cueType == nCue);

        % compute real psth
        [psth, raster_tensor, stderr] = preparePSTH_Manoj(ss, 'spikes', trialIndicesPerCond.cueOn.(cue_cond_field), [ss.cueOnset], timePoints, binsize_rescaled);
        PSTH.cueOn.(cue_cond_field).ss.raster_tensor = raster_tensor;
        PSTH.cueOn.(cue_cond_field).ss.psth = psth;
        PSTH.cueOn.(cue_cond_field).ss.stderr = stderr;

        % compute lfads psth
        [psth, raster_tensor, stderr] = preparePSTH_Manoj(tc, 'spikes', trialIndicesPerCond.cueOn.(cue_cond_field), [ss.cueOnset], timePoints, binsize_rescaled);
        PSTH.cueOn.(cue_cond_field).tc.raster_tensor = raster_tensor;
        PSTH.cueOn.(cue_cond_field).tc.psth = psth;
        PSTH.cueOn.(cue_cond_field).tc.stderr = stderr;
    end
end
                

%% first let's plot barOn
saveDir = fullfile(saveRoot, 'errorBar','1_5_sd_cons');
if ~isdir(saveDir)
    mkdir(saveDir)
end
cd(saveDir)
nNeurons = size(PSTH.barOn.Vert.ss.psth, 1);
%%
for n = 1:nNeurons
    f1 = figure;
    y_max = 1.2*max([PSTH.barOn.Vert.tc.psth(n, :), PSTH.barOn.Vert.ss.psth(n, :), PSTH.barOn.Hori.tc.psth(n, :), PSTH.barOn.Hori.ss.psth(n, :)]);
    s1 = subplot(3,2,1);
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(PSTH.barOn.Vert.ss.psth(n,:), 'r', 'LineWidth', 2, 'DisplayName','Vertical Bars');
    hold on
    plot(PSTH.barOn.Vert.ss.psth(n,:)+PSTH.barOn.Vert.ss.stderr(n,:), 'r', 'LineWidth', 1);
    hold on
    plot(PSTH.barOn.Vert.ss.psth(n,:)-PSTH.barOn.Vert.ss.stderr(n,:), 'r', 'LineWidth', 1);
    hold on
    plot(PSTH.barOn.Hori.ss.psth(n,:), 'b', 'LineWidth', 2, 'DisplayName','Horizontal Bars');
    hold on
    plot(PSTH.barOn.Hori.ss.psth(n,:) + PSTH.barOn.Hori.ss.stderr(n,:), 'b', 'LineWidth', 1);
    hold on
    plot(PSTH.barOn.Hori.ss.psth(n,:) - PSTH.barOn.Hori.ss.stderr(n,:), 'b', 'LineWidth', 1);
    hold on
    set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
    set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
    title(s1, 'SS', 'FontSize', 16);
    ylabel('Firing Rate');
    axes(s1)
    ylim([0 y_max])
    xlim([0 801]) %need to change there
    set(gca, 'fontsize', 12);    

    s2 = subplot(3,2,2);
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(PSTH.barOn.Vert.tc.psth(n,:), 'r', 'LineWidth', 2, 'DisplayName','Vertical Bars');
    hold on
    plot(PSTH.barOn.Vert.tc.psth(n,:)+PSTH.barOn.Vert.tc.stderr(n,:), 'r', 'LineWidth', 1);
    hold on
    plot(PSTH.barOn.Vert.tc.psth(n,:)-PSTH.barOn.Vert.tc.stderr(n,:), 'r', 'LineWidth', 1);
    hold on
    plot(PSTH.barOn.Hori.tc.psth(n,:), 'b', 'LineWidth', 2, 'DisplayName','Horizontal Bars');
    hold on
    plot(PSTH.barOn.Hori.tc.psth(n,:) + PSTH.barOn.Hori.tc.stderr(n,:), 'b', 'LineWidth', 1);
    hold on
    plot(PSTH.barOn.Hori.tc.psth(n,:) - PSTH.barOn.Hori.tc.stderr(n,:), 'b', 'LineWidth', 1);
    hold on
    set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
    set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
    title(s2, 'TC', 'FontSize', 16);
    ylabel('Firing Rate');
    axes(s2)
    ylim([0 y_max])
    xlim([0 801]) % need to change here to not hard coded
    set(gca, 'fontsize', 12);    
    %legend('show')
    %set(legend, 'Location', 'southeast');

    s3 = subplot(3,2,3);
    a = squeeze(PSTH.barOn.Vert.ss.raster_tensor(n,:,:));
    imagesc(a)
    set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
    set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
    set(gca, 'fontsize', 12);    
    title(s3, 'Vert');
    ylabel('Trials');

    s4 = subplot(3,2,4);
    a = squeeze(PSTH.barOn.Vert.tc.raster_tensor(n,:,:));
    imagesc(a)
    set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
    set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
    set(gca, 'fontsize', 12);    
    title(s4, 'Vert');
    ylabel('Trials');

    s5 = subplot(3,2,5);
    a = squeeze(PSTH.barOn.Hori.ss.raster_tensor(n,:,:));
    imagesc(a)
    set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
    set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
    set(gca, 'fontsize', 12);    
    title(s5, 'Hori');
    ylabel('Trials');

    s6 = subplot(3,2,6);
    a = squeeze(PSTH.barOn.Hori.tc.raster_tensor(n,:,:));
    imagesc(a)
    set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
    set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
    set(gca, 'fontsize', 12);    
    title(s6, 'Hori');
    ylabel('Trials');
    
    suptitle(['Multi-unit ' int2str(n)]);
    set(f1, 'Position', [428 200 1215 811]);
    print(f1,['Multi-unit ' int2str(n)], '-dpng');
    %printpdf(f1,int2str(nIndices(n)) )
    close;
end




