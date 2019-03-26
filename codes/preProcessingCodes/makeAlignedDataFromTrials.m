%% add your paths here.

% add paths for Feng
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/postAnalysisCodes')

%% set up names for datasets
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

%% load pre-processed trialized data (made from continuous dataset))
loadpath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
            'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/'];

% iterate over days, load each day and add it to a 'olapChopped' cell array
clear olapChopped
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_cueOnArrayOnTargetDim_HoldRel.mat', loadpath, datasets( nday ).shortName );

    tmp = load( fname );
    olapChopped{ nday } = tmp.combinedData;
end

%% load the UE data (ryan's stuff)

for nday = 1: numel( datasets )
    disp( sprintf( 'loading UEs day %g / %g', nday, numel( datasets ) ) );
    loaddir = sprintf('/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v12/M%s/MUA_GRATINGS/', datasets( nday ).longName );
    searchPattern = sprintf( '%sM%s*-evokedSpiking-*.mat', loaddir, datasets( nday ).longName );
    tmp = dir( searchPattern );
    disp( tmp(1).name );

    fname = sprintf( '%s%s', loaddir, tmp(1).name );
    data = load(fname);
    UEs{nday} = data.UE;
end

%%
rawSampleRate = 1000;
cueTimePoints = -300:900;
arrayTimePoints_hold = -400:800;
arrayTimePoints_release = -500:700;
dimTimePoints = -600:600;

%% load the UE data (ryan's stuff)
alf = {};
for nday = 1: numel( datasets )
    
    % get all the timing info from the UE struct on a per-trial basis
    %% get cueOn, arrayOn, targetDim timing from real data
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
    hold_ind = 0;
    trialNum = numel( cueStart );
    %    holdTrialNum = numel( dimStart );
    cueLocs = unique(UEs{nday}.cueLoc);
    for ntr = 1:numel( cueInds )
        % first load in cueOnset trials
        alf{ nday }( ntr ).spikeCounts = olapChopped{ nday }.r.r( ntr ).spikes( :, cueStart( ntr ) + cueTimePoints );
        alf{ nday }( ntr ).externalInputs = olapChopped{ nday }.r.r( ntr ).externalInputs( :, cueStart( ntr ) + cueTimePoints );        
        alf{ nday }( ntr ).condition = find(cueLocs == UEs{ nday }.cueLoc( ntr ));

        % second load in arrayOnset trials
        if UEs{ nday }.isHoldTrial(ntr)
            alf{ nday }( ntr + trialNum ).spikeCounts = olapChopped{ nday }.r.r( ntr ).spikes( :, arrayStart( ntr ) + arrayTimePoints_hold );
            alf{ nday }( ntr + trialNum ).externalInputs = olapChopped{ nday }.r.r( ntr ).externalInputs( :, arrayStart( ntr ) + arrayTimePoints_hold );        
            alf{ nday }( ntr + trialNum ).condition = find(cueLocs == UEs{ nday }.cueLoc( ntr )) + numel( cueLocs ); % hold trials conditions: 3 and 4
        else
            alf{ nday }( ntr + trialNum ).spikeCounts = olapChopped{ nday }.r.r( ntr ).spikes( :, arrayStart( ntr ) + arrayTimePoints_release );
            alf{ nday }( ntr + trialNum ).externalInputs = olapChopped{ nday }.r.r( ntr ).externalInputs( :, arrayStart( ntr ) + arrayTimePoints_release );        
            alf{ nday }( ntr + trialNum ).condition = find(cueLocs == UEs{ nday }.cueLoc( ntr )) + 2*numel( cueLocs ); % release trials conditions: 5 and 6
        end

        % third load in targetDim trials
        if UEs{ nday }.isHoldTrial( ntr )
            hold_ind = hold_ind + 1;
            alf{ nday }( 2*trialNum + hold_ind ).spikeCounts = olapChopped{ nday }.r.r( ntr ).spikes( :, dimStartAll( ntr ) + dimTimePoints );
            alf{ nday }( 2*trialNum + hold_ind ).externalInputs = olapChopped{ nday }.r.r( ntr ).externalInputs( :, dimStartAll( ntr ) + dimTimePoints );        
            alf{ nday }( 2*trialNum + hold_ind ).condition = find(cueLocs == UEs{ nday }.cueLoc( ntr )) + 3*numel( cueLocs ); % hold trials conditions: 7 and 8
        end
        
        
        %        alf{ nday }( ntr ).cueOnset =  cueStart( ntr );
        %        alf{ nday }( ntr ).arrayOnset = arrayStart( ntr );
        %        alf{ nday }( ntr ).arrayDim =  dimStartAll( ntr );
    end
end

%%
saveDir = '/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/preAligned/multi-day_CoAoTdHoldRel_JanToApr/withExternalInputs_twoLoc_181023';
cd(saveDir);
for nday = 1: numel(datasets)
    R = alf{nday};
    saveName = [datasets( nday ).shortName, '_cueOnArrayOnTargetDim_HoldRel.mat'];
    save(saveName, 'R');
    clear R;
end

%%