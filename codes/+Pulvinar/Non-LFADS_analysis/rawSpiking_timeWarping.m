% -------------------
% data preprocessing
% -------------------
%% add appropriate paths
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')
%addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/postAnalysisCodes')

%% name datasets
datasets(1).shortName = '170127';
datasets(1).longName = '20170127';
datasets(2).shortName = '170130';
datasets(2).longName = '20170130';
datasets(3).shortName = '170201';
datasets(3).longName = '20170201';
datasets(4).shortName = '170211';
datasets(4).longName = '20170211';
datasets(5).shortName = '170308';
datasets(5).longName = '20170308';
datasets(6).shortName = '170311';
datasets(6).longName = '20170311';

%% % load all the chopped and recombined data (post-LFADS)
loadpath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
            'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/'];
loadpath_rmCorr = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
            'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/180614data_rm_highCorr/'];
% iterate over days, load each day and add it to a 'olapChopped' cell array
clear olapChopped
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_cueOnArrayOnTargetDim_HoldRel.mat', loadpath_rmCorr, datasets( nday ).shortName );

    tmp = load( fname );
    olapChopped{ nday } = tmp.combinedData;
end

%% load the UE data (ryan's stuff)
for nday = 1: numel( datasets )
    disp( sprintf( 'loading UEs day %g / %g', nday, numel( datasets ) ) );
    loaddir = sprintf('/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v12/M%s/MUA_GRATINGS/', datasets( nday ).longName );
    %loaddir = sprintf('/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/M%s/Gratings/', datasets( nday ).longName );
    searchPattern = sprintf( '%sM%s*-evokedSpiking-*.mat', loaddir, datasets( nday ).longName );
    tmp = dir( searchPattern );
    disp( tmp(1).name );

    fname = sprintf( '%s%s', loaddir, tmp(1).name );
    data = load(fname);
    UEs{nday} = data.UE;
end

%% preprocessing spikes
binsize_rescaled = 8;
use_rates = false;
rawSampleRate = 1000;
trial_time_ms = 500;
trial_olap_ms = 100;

for nday = 1:numel(datasets)
        % get all real data
    trueSpikes = [olapChopped{ nday }.r.r.spikes];
    nSpikes = size(trueSpikes, 2);
    out = olapChopped{ nday }.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );
    tot_spikes = out.counts;
    % nPieces = size(tot_spikes, 1);
    %spikes = zeros(size(trueSpikes, 1), nSpikes);
    spikes(length(olapChopped{ nday }.r.r)).spikes = 1;
    % assembly
    nChops = 0;
    for r_trial = 1: length(olapChopped{ nday }.r.r)
        spikes(r_trial).spikes = zeros(size(olapChopped{ nday }.r.r(r_trial).spikes));
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
    rbinned_real = assembled_real2.binData({'spikes'}, [binsize_rescaled]);
    alf{ nday } = rbinned_real;
    
    sessStartTime = UEs{ nday }.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1)*rawSampleRate;

    % session start time - I set it to be the time when the monkey starts
    % fixation in the first trial
    sessEndTime = UEs{ nday }.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
    extraEndMs = 500;

    startInds = round(UEs{ nday }.fixationAndLeverTimes.firstEnterFixationTimesPreCue * rawSampleRate - sessStartTime);
    startInds(1) = 1;
    stopInds = round(UEs{ nday }.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice * rawSampleRate - sessStartTime);

    cueInds = round(UEs{ nday }.cueOnset * rawSampleRate - sessStartTime);
    arrayInds = round(UEs{ nday }.arrayOnset * rawSampleRate - sessStartTime);
    dimInds = round(UEs{ nday }.targetDim * rawSampleRate - sessStartTime);

    cueStart = cueInds - startInds;
    arrayStart = arrayInds - startInds;
    dimStart = dimInds - startInds(UEs{ nday }.isHoldTrial);
    % 'dimStart' is only defined for hold trials
    % make a version that's defined for all trials
    dimStartAll = nan( size( cueStart ) );
    dimStartAll( UEs{ nday }.isHoldTrial ) = dimStart;
    % add important info to the 'alf' struct

    for ntr = 1:numel( alf{ nday } )
        alf{ nday }( ntr ).cueOnset = round( cueStart( ntr ) / binsize_rescaled);
        alf{ nday }( ntr ).arrayOnset = round( arrayStart( ntr ) / binsize_rescaled);
        alf{ nday }( ntr ).arrayDim = round( dimStartAll( ntr ) / binsize_rescaled);
    end
end 
