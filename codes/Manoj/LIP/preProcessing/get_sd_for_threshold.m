%filename = '/mnt/scratch/feng/Remy_02262019_PUL_Raw.pl2';
%outfile = '/mnt/scratch/cpandar/Remy_02262019_PUL_spikeband.mat';
%
%tic;
%broadband2streamMinMax( filename, outfile )
%toc;
%% first load UEs
% load UE and get trialstruct
clear UE 
load_UE_path = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/UE_02182019.mat';
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

%%
spikeBandFile = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/Remy_02182019_LIP_spikeband.mat';
bb = load(spikeBandFile);
bb = bb.spikeband;
% remove common average
bb.minSpikeBand_rmCA = bb.minSpikeBand - mean(bb.minSpikeBand,2);

%% get var and remove NaN, also remove NaN for minSpikeBand
for ich = 1:size(bb.minSpikeBand,2)
    chVar{ ich } = bb.meanSquared( bb.meanSquaredChannel == ich );
    whereNan = find( isnan( chVar{ ich } ) );
    chVar{ ich } = chVar{ ich }( 1 : (whereNan(1) - 1) );
    whereNan_msb = find( isnan( bb.minSpikeBand_rmCA( : , ich ) ) );
    chMsb{ ich } = bb.minSpikeBand_rmCA(1:(whereNan_msb(1) - 1), ich);
end

%% compute changing std
chStdVec = zeros(1,32);
for ich = 1:numel( chVar )
    % % use original variance (changing value across a session)
    chStd{ich} = sqrt(chVar{ich});
    chStd{ich} = repelem(chStd{ich}, 32);
    if numel(chStd{ich}) < numel(chMsb{ich})
        num_diff = numel(chMsb{ich}) - numel(chStd{ich});
        elemToAdd = chStd{ich}(end,1)*ones(num_diff, 1);
        chStd{ich} = [chStd{ich};elemToAdd];
    else
        chStd{ich} = chStd{ich}(1:length(chMsb{ich}));
    end 
end

%% compute constant std
chStdVec = zeros(1,32);
for ich = 1:numel( chVar )
    %use mean of variance across entire session for each channel
    chStd_cons{ich} = sqrt( mean( chVar{ich} ) );
    chStdVec(ich) = chStd_cons{ich};
end

%% remove common average
bb.minSpikeBand_rmCA = bb.minSpikeBand - mean(bb.minSpikeBand,2);

%% start iteration on the # of sd to get spikes
sd_factors = [1.25, 1.5, 1.75, 2, 2.25];
for i = 1:numel(sd_factors)
    sd_factor = sd_factors(i);
    spikes = sparse(size(chMsb{1}, 1), numel(chMsb));
    for ich = 1:size(bb.minSpikeBand,2)
        chMean = mean(chMsb{ich});
        chThres = chMean - sd_factor*chStd_cons{ich};
        leftValue = chMsb{ich} - chThres;
        spikes(leftValue <= 0, ich) = 1;
    end
    % remove spikes if show on more than 25 TCs (~80%)
    allSpikesMS = sum(spikes,2);
    spikes((allSpikesMS > 25), :) = 0;
    % % get rid of MU 1 - 4
    %spikes = spikes(:, 5:end);
    stream.spikes = sparse(spikes);

    % preprocess spikes
    clear C r tc_r tc_rebinned tc
    dtMS = 1;
    C = Continuous.Continuous(stream, dtMS);
    sigma_neural = 10;
    C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );
    r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );
    tc_r = R.Rstruct(r.r);
    %
    sigma = 50;
    binsize_rescaled = 10;
    %
    tc_r.smoothFieldInR('spikes', 'spike_smoothed', sigma, 1);
    tc_rebinned = tc_r.binData({'spike_smoothed'}, [binsize_rescaled]);

    barStart = UE.barOn - UE.fixOn;
    cueStart = UE.cueOn - UE.fixOn;
    targetStart = UE.targetOn - UE.fixOn;
    tc = [];
    for ntr = 1:numel(tc_rebinned)
        tc(ntr).spikes = tc_rebinned(ntr).spike_smoothed;
        tc( ntr ).barOnset = round( barStart( ntr ) /binsize_rescaled);
        tc( ntr ).cueOnset = round( cueStart( ntr ) / binsize_rescaled);
        tc( ntr ).targetStart = round( targetStart( ntr ) / binsize_rescaled);
    end

    % get PSTHs and rasters
    events = {'barOn', 'cueOn', 'targetOn'};
    conditions = struct;
    conditions.barOn = {'Vert', 'Hori'};
    conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
    window = round([-400  400] / binsize_rescaled);
    timePoints = window(1):window(2);

    %
    for nBar = 1:numel(conditions.barOn)
        trialIndicesPerCond.barOn.(conditions.barOn{nBar}) = UE.barType == nBar;
        % compute tc psth
        [psth, raster_tensor, stderr] = preparePSTH_Manoj(tc, 'spikes', trialIndicesPerCond.barOn.(conditions.barOn{nBar}), [tc.barOnset], timePoints, binsize_rescaled);
        PSTH(i).barOn.(conditions.barOn{nBar}).tc.raster_tensor = raster_tensor;
        PSTH(i).barOn.(conditions.barOn{nBar}).tc.psth = psth;
        PSTH(i).barOn.(conditions.barOn{nBar}).tc.stderr = stderr;
        
        for nCue = 1 : numel(conditions.cueOn)
            cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
            trialIndicesPerCond.cueOn.(cue_cond_field) = (UE.barType == nBar) & (UE.cueType == nCue);

            % compute lfads psth
            [psth, raster_tensor, stderr] = preparePSTH_Manoj(tc, 'spikes', trialIndicesPerCond.cueOn.(cue_cond_field), [tc.cueOnset], timePoints, binsize_rescaled);
            PSTH(i).cueOn.(cue_cond_field).tc.raster_tensor = raster_tensor;
            PSTH(i).cueOn.(cue_cond_field).tc.psth = psth;
            PSTH(i).cueOn.(cue_cond_field).tc.stderr = stderr;
        end
    end
end

%% find top 5 channels that have highest FR
[~, sidx] = sort(mean(PSTH(3).barOn.Vert.tc.psth, 2), 'descend');
%nn = [1:5, 8:4:32];
nn = [25:32]
N = sidx(nn)'; % neuron
%N = sidx(1:10)'; % neuron
pre_post_times = [400, 400] / binsize_rescaled;
t = pre_post_times;
t = [-t(1):-1 0:t(2)] * binsize_rescaled; %t(end)=[];

%%
savePSTHdir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/PSTHs/selecting_threshold/02182019/tmp_comAve_rm_2'
if ~isdir(savePSTHdir)
    mkdir(savePSTHdir)
end
cd(savePSTHdir)
whichEvent = 'cueOn';
cond_1 = 'Exo_TL_Vert';
cond_2 = 'Exo_BL_Vert';
numS = 5;
figure
for n = N
    clf;
    i=0;
    for s = 1:5
        p_1 = PSTH(s).(whichEvent).(cond_1).tc.psth(n, :);
        p_2 = PSTH(s).(whichEvent).(cond_2).tc.psth(n, :);
        r_1 = squeeze(PSTH(s).(whichEvent).(cond_1).tc.raster_tensor(n, :, :));
        r_2 = squeeze(PSTH(s).(whichEvent).(cond_2).tc.raster_tensor(n, :, :));
        %        keyboard

        subplot(3,numS,i+1)
        plot(t, p_1, '-b')
        hold on
        plot(t, p_2, '-r')
        axis tight
        if s == 1
            yl = ylim;
            yl = [0 yl(2)*1.5];
        end
        plot([0  0], yl, '--k')
        if s == 1
            legend(cond_1, cond_2, whichEvent)
            legend('Location', 'best')
        end
        ylim(yl)
        title(sprintf('Neuron: %d - Num of sd: %d', n, sd_factors(s)));
        subplot(3,numS,i+numS+1)
        imagesc(t,[], r_1)
        if s==1
            cl1=caxis;
        end
        caxis(cl1)
        hold on
        plot([0 0], ylim, '--w')
        ylabel(cond_1)
        subplot(3,numS,i+2*numS+1)
        imagesc(t,[],r_2)
        if s==1
            cl2=caxis;
        end
        caxis(cl2)
        hold on
        plot([0 0], ylim, '--w')
        ylabel(cond_2)
        xlabel('(ms)')
        i = i + 1;
    end
    set(gcf, 'Position', [100, 10, 1800, 800])
    tightfig(gcf);
    print(gcf,['Neuron ' int2str(n)], '-dpng');
    %    print(gcf, fullfile(save_path, sprintf('PSTH_N%d_%s', n, fig_tag)), '-dpng');
end


