%% set your paths
addpath('/snel/home/cpandar/analysis_tools/');

%% 
ddir = '/snel/share/share/data/Miller_Perich/';
matname = 'Chewie_20161011_CO_FF_BL_001_stripped.mat';
%matname = 'Chewie_20150319_CO_CS_BL_001_stripped.mat';

fname = fullfile( ddir, matname );

%%
[C, startInds, stopInds, trialstruct] = Perich_loadAndProcess( fname );

%%
% smooth the spiketrains with a gaussian kernel
C.smoothField( 'spikes', 'y_smoothed', 100 );

%%
clear r R
% turn into a trialized (R) struct
r = R.Rstruct( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

% trim to successful trials only
outcomes = [r.r.result];
r.r = r.r( outcomes == 'R' );

%%
for channelToTest = 1:242
    rcopy = r.copy();
    
    % assume our "factors" arethe smoothed neural data of neurons 1-241
    % fit a GLM for neuron 242

    keepChannels = true(242, 1);
    keepChannels( channelToTest ) = false;
    keepChannels = find( keepChannels );

    for itrial = 1:numel( rcopy.r )
        rcopy.r( itrial ).y_smoothed = rcopy.r( itrial ).y_smoothed( keepChannels, : );
        rcopy.r( itrial ).z( 1, : ) = rcopy.r( itrial ).spikes( channelToTest, : );
    end

    %%
    %  bin the data
    rbinned = rcopy.binData( { 'z', 'y_smoothed' }, 200 );

    %%
    %  fit a GLM model
    numTrials = numel(rbinned);
    trainTrials = true( size( rbinned ) );
    validTrials = 1:5:numTrials;
    trainTrials( validTrials ) = false;
    trainTrials = find( trainTrials );

    model = GLM.fitGLM( rbinned( trainTrials ), 'y_smoothed', 'z', 1 );

    %%
    % test our GLM model
    dataOut = GLM.evalGLM( model, rbinned( validTrials ), 'y_smoothed', 'z', 1, 'firingRate' );

    %%
    % plot the predicted firing rates vs the actual binned spiking activity
    binnedSpikes = [ dataOut.z ];
    predictedFR = [ dataOut.firingRate ];

    h = scatter( predictedFR, binnedSpikes, '.' );
    set(h,'markeredgealpha', 0.25 );

    title( sprintf( 'Neuron %g', channelToTest ) );
    keyboard
end