%% % datasets
% Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');

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



%% % load all the prealigned data (post-LFADS)
loadpath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/preAligned/multi-day_CoAoTdHoldRel_JanToApr/withExternalInputs_twoLoc_181023/'];

% iterate over days, load each day and add it to a 'olapChopped' cell array
clear preAligned
for nday = 1:numel( datasets )
    disp( sprintf( 'loading prealigned day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_cueOnArrayOnTargetDim_HoldRel.mat', loadpath, datasets( nday ).shortName );
    tmp = load( fname );
    for itrial = 1 : numel( tmp.R )
        tmp.R(itrial).spikeCounts = full(tmp.R(itrial).spikeCounts);
        tmp.R(itrial).externalInputs = full(tmp.R(itrial).externalInputs);
    end
    alf{ nday } = tmp.R;
end

%% load rf info
loadpath_rf = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
            'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/'];

% iterate over days, load each day and add it to a 'olapChopped' cell array
clear rfLoc
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_cueOnArrayOnTargetDim_HoldRel.mat', loadpath_rf, datasets( nday ).shortName );
    tmp = load( fname );
    rfLoc{ nday } = tmp.combinedData.r.r(1).rfloc;
end

%%
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

%% load LFADS posterior means for all days
r_id = 1;
run = rc2.runs(r_id);
run.loadSequenceData(); % load sequence data in that run
run.loadPosteriorMeans(); % load posterior mean in that run
run.addPosteriorMeansToSeq();

%%

for nday = 1 : numel( datasets )
    thisDayReal = R.Rstruct( alf{ nday } );
    % first get binned spiking
    rbinned_real = thisDayReal.binData({'spikeCounts'}, [binsize_rescaled]);
    % second get binned smoothed spiking
    thisDayReal.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
    rbinned_smoothed = thisDayReal.binData({'spike_smoothed'}, [binsize_rescaled]);

    for ntr = 1:numel( alf{ nday } )
        alf{ nday }( ntr ).binned_spikes = rbinned_real( ntr ).spikeCounts;
        alf{ nday }( ntr ).binSmoothed_spikes = rbinned_smoothed( ntr ).spike_smoothed;
        
        rdim = size( run.sequenceData{ nday }(ntr).factors );
        tmp1 = run.sequenceData{ nday }(ntr).factors(:, 1:floor( rdim( 2 ) / rfactor ) * rfactor );
        tmp1 = reshape( tmp1, rdim(1), rfactor, floor( rdim( 2 ) / rfactor ) );
        tmp1 = squeeze( mean( tmp1, 2 ) );
        alf{ nday }(ntr).factors = tmp1;

        rdim2 = size( run.sequenceData{ nday }( ntr ).rates );
        tmp2 = run.sequenceData{ nday }( ntr ).rates(:, 1:floor( rdim2( 2 ) / rfactor ) * rfactor );
        tmp2 = reshape( tmp2, rdim2(1), rfactor, floor( rdim2( 2 ) / rfactor ) );
        tmp2 = squeeze( mean( tmp2, 2 ) );
        alf{ nday }(ntr).rates = tmp2;        
    end
end