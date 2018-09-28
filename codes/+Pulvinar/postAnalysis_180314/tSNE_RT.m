%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/preAligned/' ...
    'multi-day_CoAoTdHoldRel_JanToApr/withGoodNeurons_HoldRelSepForAO'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRelSepForAO_JanToApr';

% add individual datasets
Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170320_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170324_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170327_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170329_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170331_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'withGoodNeurons_Run_20180314', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180314;

% this script defines the run params
Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
for nrun = 1:numel( par6 )
    rc2.addParams( par6( nrun ) );
end
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

return;

%% Post Analysis

%% Post-running analysis - loading data and the output of LFADS
for nData = 1:length(dc.datasets)
    realData = dc.datasets(nData).loadData();
    r_real(nData) = R.Rstruct(realData.R);
end
% r_real = dc.datasets(1).loadData(); % get the original dataset (for all neurons)
% r_real = R.Rstruct(r_real.R); % put the dataset into R struct class




%% loading data and the output of LFADS
%for r_id = 1:length(rc2.runs)
for r_id = 2
    run = rc2.runs(r_id); % pull out run information
    run.loadSequenceData(); % load sequence data in that run
    run.loadPosteriorMeans(); % load posterior mean in that run
    run.addPosteriorMeansToSeq();
%     for nData = 1:length(run.sequenceData)
    nData = 11;
    RunID(r_id).r_lfads(nData) = R.Rstruct(run.sequenceData{nData});
%     end
%     r_lfads(r_id) = R.Rstruct(run.sequenceData{1}); % Put sequence data into a struct
end

%% Select the run and day you want to analyse
r_realCopy = r_real(11).copy();
r_lfadsWhole = RunID(2).r_lfads(11).copy();

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(r_realCopy.r); % get trial number
nTimesRaw = size(r_realCopy.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(r_realCopy.r(1).spikeCounts, 1); % get neuron nubmer
nTimesLFADS = size(r_lfadsWhole.r(1).rates,2);% get trial length for rebinned data that was operated by LFADS 
% modify this line if nTimes for different trials or runs are different.
nFactors = size(r_lfadsWhole.r(1).factors, 1);


%% pull out trials that will be analyzed
r_lfadsCopy = r_lfadsWhole.copy();
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


%% t-SNE analysis on response time, with raw spiking data - coarse separation on RT
savedirFive = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'postAnalysis/withGoodNeurons_Run_20180314/tSNE_onRT/tSNE_smoothedRealRate/fastVsSlow/170331/'];

timeLagAll = -20:2:20
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
    f1 = figure
    scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
    hold on
    scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
    ylabel('2nd Dimension');
    xlabel('1st Dimension');
    legend('show')
    title(['TargetDim ' num2str(timeLag*par.spikeBinMs, '%+0.0f') 'ms']);
    print(f1,['TargetDim ' num2str(timeLag*par.spikeBinMs, '%+0.0f') 'ms'], '-dpng')
    close
    
    
end


%% t-SNE analysis on response time, with LFADS factors - coarse separation on RT
savedirFive = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'postAnalysis/withGoodNeurons_Run_20180314/tSNE_onRT/tSNE_LFADSFactors/fastVsSlow/170329/'];


timeLagAll = -20:2:20
for nLag = 1:length(timeLagAll)


    timeLag = timeLagAll(nLag);
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = TimeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nFactors);
    rtVector = [realAlign.rt];
    % put all the reaction time into a vector
    clVector = [realAlign.cueLoc]; 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = lfadsAlign(i).factors(:,TimeForAnalysis);
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
    f1 = figure
    scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
    hold on
    scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
    ylabel('2nd Dimension');
    xlabel('1st Dimension');
    legend('show')
    title(['TargetDim ' num2str(timeLag*par.spikeBinMs, '%+0.0f') 'ms']);
    print(f1,['TargetDim ' num2str(timeLag*par.spikeBinMs, '%+0.0f') 'ms'], '-dpng')
    close
    
    
end
%%




%% t-SNE analysis on response time, with LFADS rates - coarse separation on RT
savedirFive = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'postAnalysis/withGoodNeurons_Run_20180314/tSNE_onRT/tSNE_LFADSRates/fastVsSlow/170329/'];


timeLagAll = -20:2:20
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
    f1 = figure
    scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
    hold on
    scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
    ylabel('2nd Dimension');
    xlabel('1st Dimension');
    legend('show')
    title(['TargetDim ' num2str(timeLag*par.spikeBinMs, '%+0.0f') 'ms']);
    print(f1,['TargetDim ' num2str(timeLag*par.spikeBinMs, '%+0.0f') 'ms'], '-dpng')
    close
    
    
end
%%














% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_TargetDimRT'];

r_lfadsCopy = r_lfads(1).copy();
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