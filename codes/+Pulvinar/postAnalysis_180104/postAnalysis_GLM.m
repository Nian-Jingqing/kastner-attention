%% build the dataset collection
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
datasetPath = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/datasets'];



%% Locate and specify the datasets
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'CO_AO_TD_HoldRel20170608';

% add individual datasets
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel_001.mat');
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel_002.mat');
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel_003.mat');
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel_004.mat');
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel_005.mat');
Pulvinar.Dataset(dc, 'cueOnArrayOnTargetDim_HoldRel.mat');
% MyExperiment.Dataset(dc, 'dataset002.mat');%Example code for lorenz
% example - for future use
% MyExperiment.Dataset(dc, 'dataset003.mat');

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = '/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/CueOnArrayOnTargetDim_HoldRel/runs';
rc = Pulvinar.RunCollection(runRoot, 'secondRun_20180104', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc.version = 20180104;

%% Set parameters for the entire run collection

par = Pulvinar.RunParams;
par.spikeBinMs = 4; % rebin the data at 4 ms
par.c_co_dim = 3; % number of units in controller --> 3-d inputs to generator
par.c_batch_size = 210; % must be < 1/5 of the min trial count % total trial number is 1101, so I chose this to be 210
par.c_factors_dim = 20; % and manually set it for multisession stitched models
par.useAlignmentMatrix = false; % use alignment matrices initial guess for multisession stitching

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_ci_enc_dim = 64; % number of units in encoder for controller input
par.c_con_dim = 64; % number of units in controller

par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo

% add a single set of parameters to this run collection. Additional
% parameters can be added. LFADS.RunParams is a value class, unlike the other objects
% which are handle classes, so you can modify par freely.
rc.addParams(par);

%% Add RunSpecs

% Run a single model for each dataset, and one stitched run with all datasets

% add each individual run
for iR = 1:dc.nDatasets
    runSpec = Pulvinar.RunSpec(dc.datasets(iR).getSingleRunName(), dc, dc.datasets(iR).name);
    rc.addRunSpec(runSpec);
end

% add the final stitching run with all datasets
%rc.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

% adding a return here allows you to call this script to recreate all of
% the objects here for subsequent analysis after the actual LFADS models
% have been trained. The code below will setup the LFADS runs in the first
% place.

return;



%% Post-running analysis - loading data and the output of LFADS
r_real = dc.datasets(end).loadData(); % get the original dataset (for all neurons)
r_real = R.Rstruct(r_real.R); % put the dataset into R struct class
for r_id = 1:length(dc.datasets)
    run = rc.runs(r_id); % pull out run information
    run.loadSequenceData(); % load sequence data in that run
    run.loadPosteriorMeans(); % load posterior mean in that run
    run.addPosteriorMeansToSeq();
    r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct
end

%% get experiment info (nTrials, nTimes, nNeurons)
nTrials = length(r_real.r); % get trial number
nTimesRaw = size(r_real.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_real.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfads(1).r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfads(1).r(6).factors, 1);

%% GLM analysis
savedirOne = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/GLM_analysis/Rebin_200ms'];
rng('default');
for r_id = 1:5
    saveFolders = dir(savedirOne);
    savedirOne_sub = fullfile(savedirOne, saveFolders(r_id+2).name);  
    cd(savedirOne_sub);
    rGLM = r_real.copy();
    x = randperm(nNeurons);
    keepChannels = x(1:75);
    channelToTest = x(76:end);
    for itrial = 1:numel( rGLM.r )
        rGLM.r( itrial ).factors = r_lfads(r_id).r( itrial ).factors;
        rGLM.r( itrial ).heldout = rGLM.r( itrial ).spikeCounts( channelToTest, : );
    end
    
    %% bin the data
    rbinned = rGLM.binData({'heldout', 'factors'}, [200, 50]);
    % bin size 20 for 'heldout' spiking and 5 for 'factor'
    % bin size = 20 ms.
    
    %% fit a GLM model
    trainTrials = true( size( rbinned ) );
    validTrials = 1:5:nTrials;
    trainTrials( validTrials ) = false;
    trainTrials = find( trainTrials );
    predictedFR = zeros( numel(channelToTest), numel(validTrials)*size(rbinned(1).heldout,2));
    % initialize a matrix to store predicted FR for heldout neurons.
    % size: 31 x (221*80)
    for ineuron = 1:numel(channelToTest)
        model = GLM.fitGLM( rbinned( trainTrials ), 'factors', 'heldout', ineuron );
        dataOut = GLM.evalGLM( model, rbinned( validTrials ), 'factors', 'heldout', ineuron, 'firingRate' );
        predictedFR(ineuron,:) = [ dataOut.firingRate ];
    end
    binnedSpikes = [dataOut.heldout];
    
    %% Plot the predicted firing rates vs the actual binned spiking activity
    for ineuron = 1:numel(channelToTest)
        f1 = figure;
        scatter( predictedFR(ineuron,:), binnedSpikes(ineuron,:), 10, 'b','filled');
%         set(gca,'markeredgealpha', 0.25 );
        title( sprintf( 'Neuron %g', ineuron ) );
        ylabel('actual spiking');
        xlabel('predicted FR');
        print(f1,['Neuron ' int2str(ineuron)], '-dpng');
        close;
    end
    
%     %% shuffle test
%     s = randperm(31);
%     shuffled_spiking = binnedSpikes(s,:);
%     for ineuron = 1:numel(channelToTest)
%         f2 = figure;
%         scatter( predictedFR(ineuron,:), shuffled_spiking(ineuron,:), 10, 'b','filled');
% %         set(gca,'markeredgealpha', 0.25 );
%         title( sprintf( 'Shuffled Neuron %g', ineuron ) );
%         ylabel('actual spiking');
%         xlabel('predicted FR');
%         print(f2,['Neuron ' int2str(ineuron)], '-dpng');
%         close;
%     end
end

%% loading and put real data into Rstruct
r_real = dc.datasets(end).loadData(); % get the original dataset (for all neurons)
r_real = R.Rstruct(r_real.R); % put the dataset into R struct class

%% get dataset info

nTrials = length(r_real.r); % get trial number
nTimesRaw = size(r_real.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_real.r(1).spikeCounts, 1); % get neuron nubmer
%% smooth the raw spiking
sigma = 100;
r_real.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);

%%
savedirOne = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/GLM_analysis/HeldoutSpikingVersusSmoothedFiring/Rebin_200ms'];
cd(savedirOne);
for channelToTest = 1:106
    rcopy = r_real.copy();
    
    % assume our "factors" are the smoothed neural data of neurons 1-105
    % fit a GLM for neuron 106

    keepChannels = true(106, 1);
    keepChannels( channelToTest ) = false;
    keepChannels = find( keepChannels );

    for itrial = 1:numel( rcopy.r )
        rcopy.r( itrial ).spike_smoothed = rcopy.r( itrial ).spike_smoothed( keepChannels, : );
        rcopy.r( itrial ).spikeCounts = rcopy.r( itrial ).spikeCounts( channelToTest, : );
    end

    %%
    %  bin the data
    rbinned = rcopy.binData( { 'spikeCounts', 'spike_smoothed' }, [200, 200] );

    %%
    %  fit a GLM model
    numTrials = numel(rbinned);
    trainTrials = true( size( rbinned ) );
    validTrials = 1:5:numTrials;
    trainTrials( validTrials ) = false;
    trainTrials = find( trainTrials );

    model = GLM.fitGLM( rbinned( trainTrials ), 'spike_smoothed', 'spikeCounts', 1 );

    %%
    % test our GLM model
    dataOut = GLM.evalGLM( model, rbinned( validTrials ), 'spike_smoothed', 'spikeCounts', 1, 'firingRate' );

    %%
    % plot the predicted firing rates vs the actual binned spiking activity
    binnedSpikes = [ dataOut.spikeCounts ];
    predictedFR = [ dataOut.firingRate ];
    f1 = figure;
    scatter( predictedFR, binnedSpikes, 10, 'b','filled' );

    title( sprintf( 'Neuron %g', channelToTest ) );
    ylabel('actual spiking');
    xlabel('predicted FR');
    print(f1,['Neuron ' int2str(channelToTest)], '-dpng');
    close;
    
end

