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

%% Post Analysis

%% %% Post-running analysis - loading data and the output of LFADS
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

%% t-SNE analysis on rate vector over all trials, color and shape by cueLoc

% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_cueOnsetCueLoc'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
alignType = 1; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType);
realAlign = r_realCopy.r(alignIx == alignType);

nTrials = length(lfadsAlign);
% nTrials is specific to alignType now

timeLagAll = -25:5:100;
for nLag = 1:length(timeLagAll)


    timeLag = timeLagAll(nLag);
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = TimeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    rtVector = [realAlign.rt];
    % put all the reaction time into a vector
    clVector = [realAlign.cueLoc]; 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = lfadsAlign(i).rates(:,TimeForAnalysis);
        % pull out the neural data for all neuron at a specific time in that trial,
        % and put them into the ith row in the rate matrix
    end
    Rate_tsne = tsne(RateMatrix);
    % Actually run tsne
    trialShape = ['o', '>', 's', 'd'];
    trialColor = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0.8 0]];
    % Use different shapes to indicate cue location of the trial.
    % circle and red for cueLoc = 1; right-pointing triangle and green for cueLoc = 2;
    % square and blue for cueLoc = 3; diamond and dark yellow for cueLoc = 4;
    cd(savedirFive)
    f5 = figure
    % plot the tsne results and add coloring based on 'trialColors' matrix
    for ntrial = 1:nTrials
        h = scatter(Rate_tsne(ntrial, 1), Rate_tsne(ntrial, 2));
        set(h, 'markerfacecolor', 'none');
        set(h, 'markeredgecolor',trialColor(clVector(ntrial),:));
        set(h, 'Marker', trialShape(clVector(ntrial)));
        % use the cued Location value in clVector to index to the proper trial
        % shape in trialShape vector
        hold on;
    end
    xlabel('1st dimension');
    ylabel('2nd dimension');
    title(['cueOnset ' num2str(timeLag*4, '%+0.0f') 'ms']);
    print(f5,['cueOnset ' num2str(timeLag*4, '%+0.0f') 'ms'], '-dpng');
    close
end


%% t-SNE analysis on real spiking over all trials, color and shape by cueLoc

% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_onNeuronFiringRateCueOnsetCueLoc'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
sigma = 50;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
rbinned = r_realCopy.binData({'spike_smoothed'}, [par.spikeBinMs]);
alignType = 1; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType);
realAlign = r_realCopy.r(alignIx == alignType);
realSmoothed = rbinned(alignIx == alignType);

nTrials = length(realAlign);
% nTrials is specific to alignType now

timeLagAll = -10:5:75;
for nLag = 1:length(timeLagAll)


    timeLag = timeLagAll(nLag);
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = TimeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    rtVector = [realAlign.rt];
    % put all the reaction time into a vector
    clVector = [realAlign.cueLoc]; 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = realSmoothed(i).spike_smoothed(:,TimeForAnalysis);
        % pull out the neural data for all neuron at a specific time in that trial,
        % and put them into the ith row in the rate matrix
    end
    Rate_tsne = tsne(RateMatrix);
    % Actually run tsne
    trialShape = ['o', '>', 's', 'd'];
    trialColor = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0.8 0]];
    % Use different shapes to indicate cue location of the trial.
    % circle and red for cueLoc = 1; right-pointing triangle and green for cueLoc = 2;
    % square and blue for cueLoc = 3; diamond and dark yellow for cueLoc = 4;
    cd(savedirFive)
    f5 = figure
    % plot the tsne results and add coloring based on 'trialColors' matrix
    for ntrial = 1:nTrials
        h = scatter(Rate_tsne(ntrial, 1), Rate_tsne(ntrial, 2));
        set(h, 'markerfacecolor', 'none');
        set(h, 'markeredgecolor',trialColor(clVector(ntrial),:));
        set(h, 'Marker', trialShape(clVector(ntrial)));
        % use the cued Location value in clVector to index to the proper trial
        % shape in trialShape vector
        hold on;
    end
    xlabel('1st dimension');
    ylabel('2nd dimension');
    title(['cueOnset ' num2str(timeLag*4, '%+0.0f') 'ms']);
    print(f5,['cueOnset ' num2str(timeLag*4, '%+0.0f') 'ms'], '-dpng');
    close
end

%% t-SNE analysis on response time
% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_TargetDimRT'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
alignType = 3; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType & holdIx == 1);
realAlign = r_realCopy.r(alignIx == alignType & holdIx == 1);

nTrials = length(lfadsAlign);
% nTrials is specific to alignType now

timeLagAll = -25:5:100;
for nLag = 1:length(timeLagAll)


    timeLag = timeLagAll(nLag);
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = TimeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    rtVector = [realAlign.rt];
    % put all the reaction time into a vector
    clVector = [realAlign.cueLoc]; 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = lfadsAlign(i).rates(:,TimeForAnalysis);
        % pull out the neural data for all neuron at a specific time in that trial,
        % and put them into the ith row in the rate matrix
    end
    Rate_tsne = tsne(RateMatrix);

%% Colormap plot based on RT
    wcolormap = colormap('winter'); % 64 colors map in a sequence from blue to green
    conversionFactor = size(wcolormap, 1) / numel(rtVector);
    % convert the number of trials to the number of color, so every trial has a
    % color
    [~, trialOrder] = sort(rtVector, 'ascend'); 
    % sort trials based on RT. Each element in "trialOrder" indicates the index
    % of that element in original rtVector. In short, this means, if you sort
    % the rtVector based on RT, trialOrder tells you that, for this RT, which
    % trial it is.
    trialColors = zeros(numel(rtVector), 3);
    % initialize a trialColors matrix to tag each trial with a specific RT
    % color
    for ntrial = 1:numel(trialOrder)
        colorIndex = ceil(ntrial * conversionFactor); 
        % get the color index for that trial (the trial sorted in the sequence from fastest to slowest RT)
        trialColors(trialOrder(ntrial), :) = wcolormap(colorIndex, :);
        % use trialOrder to index to the original trial sequence and put the
        % correct color there
    end

    cd(savedirFive);
    f5 = figure
    % plot the tsne results and add coloring based on 'trialColors' matrix
    for ntrial = 1:nTrials
        h = scatter(Rate_tsne(ntrial, 1), Rate_tsne(ntrial, 2));
        set(h, 'markerfacecolor', trialColors(ntrial, :));
        set(h, 'markeredgecolor','none');
%         set(h, 'Marker', trialShape(clVector(ntrial)));
        % use the cued Location value in clVector to index to the proper trial
        % shape in trialShape vector
        hold on;
    end
    xlabel('1st dimension');
    ylabel('2nd dimension');
    title(['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms']);
    print(f5,['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms'], '-dpng');
    close
end

%% t-SNE analysis on neuron firing rate - reaction time

% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_onNeuronFiringRateTargetDimRT'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
sigma = 50;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
rbinned = r_realCopy.binData({'spike_smoothed'}, [par.spikeBinMs]);
alignType = 3; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType & holdIx == 1);
realAlign = r_realCopy.r(alignIx == alignType & holdIx == 1);
realSmoothed = rbinned(alignIx == alignType & holdIx == 1);


nTrials = length(realAlign);
% nTrials is specific to alignType now

timeLagAll = -50:5:50;
for nLag = 1:length(timeLagAll)


    timeLag = timeLagAll(nLag);
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = TimeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    rtVector = [realAlign.rt];
    % put all the reaction time into a vector
    clVector = [realAlign.cueLoc]; 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = realSmoothed(i).spike_smoothed(:,TimeForAnalysis);
        % pull out the neural data for all neuron at a specific time in that trial,
        % and put them into the ith row in the rate matrix
    end
    Rate_tsne = tsne(RateMatrix);

%% Colormap plot based on RT
    wcolormap = colormap('winter'); % 64 colors map in a sequence from blue to green
    conversionFactor = size(wcolormap, 1) / numel(rtVector);
    % convert the number of trials to the number of color, so every trial has a
    % color
    [~, trialOrder] = sort(rtVector, 'ascend'); 
    % sort trials based on RT. Each element in "trialOrder" indicates the index
    % of that element in original rtVector. In short, this means, if you sort
    % the rtVector based on RT, trialOrder tells you that, for this RT, which
    % trial it is.
    trialColors = zeros(numel(rtVector), 3);
    % initialize a trialColors matrix to tag each trial with a specific RT
    % color
    for ntrial = 1:numel(trialOrder)
        colorIndex = ceil(ntrial * conversionFactor); 
        % get the color index for that trial (the trial sorted in the sequence from fastest to slowest RT)
        trialColors(trialOrder(ntrial), :) = wcolormap(colorIndex, :);
        % use trialOrder to index to the original trial sequence and put the
        % correct color there
    end

    cd(savedirFive);
    f5 = figure
    % plot the tsne results and add coloring based on 'trialColors' matrix
    for ntrial = 1:nTrials
        h = scatter(Rate_tsne(ntrial, 1), Rate_tsne(ntrial, 2));
        set(h, 'markerfacecolor', trialColors(ntrial, :));
        set(h, 'markeredgecolor','none');
%         set(h, 'Marker', trialShape(clVector(ntrial)));
        % use the cued Location value in clVector to index to the proper trial
        % shape in trialShape vector
        hold on;
    end
    xlabel('1st dimension');
    ylabel('2nd dimension');
    title(['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms']);
    print(f5,['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms'], '-dpng');
    close
end

%% Color plot based on coarse seperation based on RT
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_TargetDimRTCoarseSeparation'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
alignType = 3; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType & holdIx == 1);
realAlign = r_realCopy.r(alignIx == alignType & holdIx == 1);

nTrials = length(lfadsAlign);
% nTrials is specific to alignType now

timeLagAll = 205:5:215;
for nLag = 1:length(timeLagAll)


    timeLag = timeLagAll(nLag);
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = TimeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    rtVector = [realAlign.rt];
    % put all the reaction time into a vector
    clVector = [realAlign.cueLoc]; 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = lfadsAlign(i).rates(:,TimeForAnalysis);
        % pull out the neural data for all neuron at a specific time in that trial,
        % and put them into the ith row in the rate matrix
    end
    Rate_tsne = tsne(RateMatrix);
    
    %% Color plot based on coarse seperation based on RT

     rtTagSlow = rtVector > median(rtVector)+0.05;
      % any rt > median rt + 0.05s in rtVector is considered to be slow, indicated by
     % logical 1
     rtTagFast = rtVector < median(rtVector)-0.05;
     % any rt < median rt - 0.05s in rtVector is considered to be fast, indicated by
     % logical 1


    rtSlowRates = RateMatrix(rtTagSlow,:);
    % Use rtTag to indice trials with slow rt
    rtFastRates = RateMatrix(rtTagFast,:);
    % Use rtTag to indice trials with fast rt
%     Rate_tsne = tsne(RateMatrix);
    % apply tsne to reduce the dimensionality of 106 neurons to 2d
    cd(savedirFive)
    f5 = figure
    scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
    hold on
    scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
    ylabel('2nd Dimension');
    xlabel('1st Dimension');
    legend('show')
    title(['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms']);
    print(f5,['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms'], '-dpng')
    close
    
    
end

%% t-SNE on neuron firing rate - Color plot based on coarse seperation based on RT
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_onNeuronFiringRateTargetDimRTCoarseSeparation'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
sigma = 50;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
rbinned = r_realCopy.binData({'spike_smoothed'}, [par.spikeBinMs]);
alignType = 3; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType & holdIx == 1);
realAlign = r_realCopy.r(alignIx == alignType & holdIx == 1);
realSmoothed = rbinned(alignIx == alignType & holdIx == 1);


nTrials = length(realAlign);

timeLagAll = -50:5:50;
for nLag = 1:length(timeLagAll)


    timeLag = timeLagAll(nLag);
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = TimeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    rtVector = [realAlign.rt];
    % put all the reaction time into a vector
    clVector = [realAlign.cueLoc]; 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = realSmoothed(i).spike_smoothed(:,TimeForAnalysis);
        % pull out the neural data for all neuron at a specific time in that trial,
        % and put them into the ith row in the rate matrix
    end
    Rate_tsne = tsne(RateMatrix);
    
    %% Color plot based on coarse seperation based on RT

     rtTagSlow = rtVector > median(rtVector)+0.05;
      % any rt > median rt + 0.05s in rtVector is considered to be slow, indicated by
     % logical 1
     rtTagFast = rtVector < median(rtVector)-0.05;
     % any rt < median rt - 0.05s in rtVector is considered to be fast, indicated by
     % logical 1


    rtSlowRates = RateMatrix(rtTagSlow,:);
    % Use rtTag to indice trials with slow rt
    rtFastRates = RateMatrix(rtTagFast,:);
    % Use rtTag to indice trials with fast rt
%     Rate_tsne = tsne(RateMatrix);
    % apply tsne to reduce the dimensionality of 106 neurons to 2d
    cd(savedirFive)
    f5 = figure
    scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
    hold on
    scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
    ylabel('2nd Dimension');
    xlabel('1st Dimension');
    legend('show')
    title(['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms']);
    print(f5,['TargetDim ' num2str(timeLag*4, '%+0.0f') 'ms'], '-dpng')
    close
    
    
end