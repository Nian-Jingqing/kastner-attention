%% build the dataset collection
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
datasetPath = ['/snel/share/share/derived/kastner/data_processed/singleSession/M20170608_PUL_all-g2-g3-g4-evokedSpiking/' ...
    'preAligned/CueOnArrayOnTargetDim_HoldRel/datasets'];



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

runRoot = ['//snel/share/share/derived/kastner/LFADS_runs/pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/runs'];
rc = Pulvinar.RunCollection(runRoot, 'secondRun_20180104', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc.version = 20180104;

Pulvinar.definePulvinarRunParams;

return;

%% Prepare LFADS input

% generate all of the data files LFADS needs to run everything
rc.prepareForLFADS();

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
%rc.writeShellScriptRunQueue('display', 50, 'maxTasksSimultaneously', 4, 'gpuList', [0 1], 'virtualenv', 'tensorflow');
rc.writeShellScriptRunQueue('display', 9, 'maxTasksSimultaneously', 4, 'gpuList', [0 1]);


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
%% things before
% rate_struct = r.sequenceData{1}; % Put sequence data into a struct

%% Post-running analysis - Plots of LFADS rates and spiking activity for example trials

savedirTwo = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/RasterVsLFADSRates/targetDim'];

cd(savedirTwo);


r_lfadsWhole = r_lfads(6).copy();
r_realCopy = r_real.copy();
sigma = 20;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
cutOff = sigma;
cutSpiking = zeros(nNeurons, nTimesRaw - 2*cutOff);
cutSmooth = zeros(nNeurons, nTimesRaw - 2*cutOff);
cutLFADSRate = zeros(nNeurons, nTimesLFADS - 2*cutOff/par.spikeBinMs);

trialSelectIx = [(1:30) (400:429) (799:828)];
r_realCopy.r = r_realCopy.r(trialSelectIx);
r_lfadsWhole.r = r_lfadsWhole.r(trialSelectIx);

for i = 61:90
    f3 = figure;
    sPlot_Spiking = subplot(3,1,1) % make the first subplot for plotting spiking activity (all neurons) for the trial
    cutSpiking = r_realCopy.r(i).spikeCounts(:, (cutOff+1):(end-cutOff));
    imagesc(cutSpiking); % "r_realCopy.r(i).spikeCounts" contains the spiking activity of all neurons for the trial
    set(gca,'XTick',[1 size(cutSpiking,2)]); 
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-780','+780'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_Spiking, 'Spiking Activity');
    ylabel('neurons');
    hold on
    plot([0.5*size(cutSpiking,2) 0.5*size(cutSpiking,2)],[1 nNeurons],'g'); % Plot a vertical line at aligned time
    
    
    sPlot_SmoothSpiking = subplot(3,1,2) % make the second subplot for plotting smoothed spiking activity (all neurons) for the trial
    cutSmooth = r_realCopy.r(i).spike_smoothed(:, (cutOff+1):(end-cutOff));
    imagesc(cutSmooth);
    set(gca,'XTick',[1 size(cutSmooth,2)]); 
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-780','+780'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_SmoothSpiking, 'Smoothed Spiking Activity');
    ylabel('neurons');
    hold on
    plot([0.5*size(cutSmooth,2) 0.5*size(cutSmooth,2)],[1 nNeurons],'g');
    
    sPlot_LFADSRates = subplot(3,1,3) % % make the third subplot for plotting LFADS rates (all neurons) for the trial
    cutLFADSRate = r_lfadsWhole.r(i).rates(:, ((cutOff/par.spikeBinMs) +1):(end-(cutOff/par.spikeBinMs)));
    imagesc(cutLFADSRate); 
    set(gca,'XTick',[1 0.5*size(cutLFADSRate,2) size(cutLFADSRate,2)]); 
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-780','CueOnset','+780'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_LFADSRates, 'LFADS Rates');
    ylabel('neurons');
    xlabel('time (ms)');
    hold on
    plot([0.5*size(cutLFADSRate,2) 0.5*size(cutLFADSRate,2)],[1 nNeurons],'g');
    
    suptitle(['Trial' int2str(i + 738)]);
    
    print(f3,['Trial ' int2str(i + 738)], '-dpng');
    close;
end


%% Post-running analysis - Compare avg neuron firing rate and avg LFADS rate for different cue locations
% prepare for plotting
r_lfadsWhole = r_lfads(6).copy();
r_realCopy = r_real.copy();
sigma = 10;
r_realCopy.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);
cutOff_smoothedSpiking = par.spikeBinMs*ceil(sigma/par.spikeBinMs);
r_realCopy.postSmoothCutOff( 'spike_smoothed', 'spike_cutOff', cutOff_smoothedSpiking);
cutOff_LFADSRates = ceil(sigma/par.spikeBinMs);
r_lfadsWhole.postSmoothCutOff( 'rates', 'rates_cutOff', cutOff_LFADSRates);
rbinned = r_realCopy.binData({'spike_cutOff'}, [par.spikeBinMs]);

AllTypeSmoothedSpikingMatrix = permute(cat(3, rbinned.spike_cutOff), [3 2 1]);
% re-organize the spiking data to nTrials x nTimes x nNeurons
AllTypeWholeLFADSRatesMatrix = permute(cat(3, r_lfadsWhole.r.rates_cutOff), [3 2 1]);
% re-organize the LFADS rates to nTrials x nTimes x nNeurons

nTimesLFADS = size(r_lfads(1).r(1).rates,2) - 2*cutOff_LFADSRates; 
% get new nTimesLFADS after cut off
AllclVector = arrayfun(@(x) x.cueLoc, r_realCopy.r); 
% pull out cue location for each trial and put into a vector
rfLoc = r_realCopy.r(1).rfloc;
% pull out receptive field for each neuron
%% separation by align types
alignType = 2; % select arrayOnset type
alignIx = [r_realCopy.r.alignType];
SmoothedSpikingMatrix = AllTypeSmoothedSpikingMatrix(alignIx == alignType,:,:);
WholeLFADSRatesMatrix = AllTypeWholeLFADSRatesMatrix(alignIx == alignType,:,:);
clVector = AllclVector(alignIx == 1);
nTrials = length(clVector);
%%

AllCueLoc = unique(clVector);
rebinSize = par.spikeBinMs; % data was rebinned from 1ms to 2ms during LFADS
AvgFiringRate = zeros(nNeurons, nTimesLFADS);
% initialize a avg firing rate matrix to store avg true firing rate (nNeurons x n Rebinned Times)
% for n = 1:nNeuron % loop over all neurons
AvgLFADSRate = zeros(nNeurons, nTimesLFADS);
% initialize a avg LFADS rate matrix to store avg LFADS rates (nNeurons x n Rebinned Times)
AvgFiringRate_InRF = zeros(nNeurons, nTimesLFADS);
AvgFiringRate_OffRF = zeros(nNeurons, nTimesLFADS);
AvgFiringRate_OutRF1 = zeros(nNeurons, nTimesLFADS);
AvgFiringRate_OutRF2 = zeros(nNeurons, nTimesLFADS);
% initialize avg Firing rates matrix for differet cue location
AvgLFADSRate_InRF = zeros(nNeurons, nTimesLFADS);
AvgLFADSRate_OffRF = zeros(nNeurons, nTimesLFADS);
AvgLFADSRate_OutRF1 = zeros(nNeurons, nTimesLFADS);
AvgLFADSRate_OutRF2 = zeros(nNeurons, nTimesLFADS);
% initialize avg LFADS rates matrix for differet cue location
TrialIndexInRF = false(nNeurons, nTrials);
TrialIndexOffRF = false(nNeurons, nTrials);
TrialIndexOutRF1 = false(nNeurons, nTrials);
TrialIndexOutRF2 = false(nNeurons, nTrials);
% initialize a matrix for each cue location to store the trial index. Ones
% would be index for trials that have the corresponding cue location

for n = 1:nNeurons % loop over all neurons

AvgFiringRate(n,:) = (sum(SmoothedSpikingMatrix(:,:,n),1))*(1/nTrials)*500;
% Store the avg firing rate to the avgFiringRate matrix
AvgLFADSRate(n,:) = (sum(WholeLFADSRatesMatrix(:,:,n),1))*(1/nTrials);

TrialIndexInRF(n,:) = clVector == rfLoc(n,1);
nInRF = length(clVector(clVector == rfLoc(n,1)));
% find the number of trials that cue loc is in RF
TrialIndexOffRF(n,:) = clVector == rfLoc(n,2);
% Find the index of the trials that cue loc is in Rf or off Rf
nOffRF = length(clVector(clVector == rfLoc(n,2)));
% find the number of trials that cue loc is Off RF
OutRF = AllCueLoc(~ismember(AllCueLoc, rfLoc(n,:))); 
% Find the cue locations that are not in or off Rf of this neuron
TrialIndexOutRF1(n,:) = clVector == OutRF(1);
% always find the index of the trials that cue loc number is smaller and
% put these index into the OutRF1
nOutRF1 = length(clVector(clVector == OutRF(1)));
% find the number of trials that cue loc is OutRF 1
TrialIndexOutRF2(n,:) = clVector == OutRF(2);
% always find the index of the trials that cue loc number is larger and
% put these index into the OutRF2
nOutRF2 = length(clVector(clVector == OutRF(2)));
% find the number of trials that cue loc is OutRF 2

AvgFiringRate_InRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF)*500;
AvgFiringRate_OffRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF)*500;
AvgFiringRate_OutRF1(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1)*500;
AvgFiringRate_OutRF2(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2)*500;

AvgLFADSRate_InRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF);
AvgLFADSRate_OffRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF);
AvgLFADSRate_OutRF1(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1);
AvgLFADSRate_OutRF2(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2);

end

%% Plots of avg neuron firing rate and compare to avg LFADS rate （no separation between cueLoc）
savedirTwo = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7' ...
    '/TargetDim/postAnalysis/AvgFiringPlusRasterNoCueLoc/Gaussian50ms'];

cd(savedirTwo);
n = 1;
for n = 1:20 % loop over all neurons
    f2 = figure
    subplot(2,1,1)
    plot(AvgLFADSRate(n,:), 'b', 'DisplayName','LFADS avg FR'); %Plot avg LFADS rate for that neuron
    hold on
    plot(AvgFiringRate(n,:), 'r', 'DisplayName','Real avg FR'); % plot avg true firing rate for that neuron
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(['Neuron' int2str(n)]);
    ylabel('Firing Rate');
    xlabel('time (ms)');
    legend('show')
%     hold on
%     TopOfTimeZero = zeros(1,2);
%     TopOfTimeZero(1,1)= max(AvgLFADSRate(n,:)); TopOfTimeZero(1,2)= max(AvgFiringRate(n,:));
%     TopOfTimeZero = max(TopOfTimeZero);
%     plot([0.5*nTimesLFADS 0.5*nTimesLFADS],[1 TopOfTimeZero],'g', 'DisplayName','Time at TargetDim'); 
%     % Plot a vertical line at aligned time
    subplot(2,1,2)
    imagesc(WholeSpikingMatrix(:,:,n));
    set(gca,'XTick',[1 0.5*nTimesRaw nTimesRaw]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title('Spiking raster');
    ylabel('Trials');
    % ylim([0 nOutRF2]);
    xlabel('time (ms)');
    
    suptitle(['Neuron' int2str(n)]);
    set(f2, 'Position', [200 100 1000 800]);
    
    
    print(f2,['Neuron ' int2str(n)], '-dpng');
    close;
end

%% plots of avg neuron firing rate and avg LFADS rate for different cue locations， plus spiking rasters
savedirThree = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/PSTH/Gaussian10ms/arrayOnset'];

cd(savedirThree);


for n = 1:20
    f2 = figure
    
    sPlot_avgRealRates = subplot(5,2,1) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
    plot(AvgFiringRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
    hold on
    plot(AvgFiringRate_OutRF2(n,:), 'y', 'DisplayName','Cue outside RF 2');
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 2
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgRealRates, 'Avg Real Firing Rate');
    ylabel('Firing Rate');
    set(gca,'xticklabel',{[]});
    
%     xlabel('time (ms)');
%     legend('show')
    
    sPlot_avgLFADSRates = subplot(5,2,2) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgLFADSRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgLFADSRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    hold on 
    plot(AvgLFADSRate_OutRF1(n,:), 'g', 'DisplayName','Cue outside RF 1');
    hold on
    plot(AvgLFADSRate_OutRF2(n,:), 'y', 'DisplayName','Cue outside RF 2');
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgLFADSRates, 'Avg LFADS Firing Rate');
    ylabel('Firing Rate');
    set(gca,'YLim',sPlot_avgRealRates.YLim);
    xlabel('time (ms)');
    Legend = legend('show');
    set(Legend,'Position',[0.7387 0.6 0.1530 0.0655]);
    
    sPlot_InRFRaster = subplot(5,2,3)
    imagesc(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'});
    set(gca,'xticklabel',{[]});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    sPlot_OffRFRaster = subplot(5,2,5)
    imagesc(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'xticklabel',{[]})
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    % ylim([0 nOffRF]);
   
    
    sPlot_OutRF1Raster = subplot(5,2,7)
    imagesc(SmoothedSpikingMatrix(TrialIndexOutRF1(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    set(gca,'xticklabel',{[]})
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OutRF1Raster, 'Cue outside RF 1');
    ylabel('Trials');
    % ylim([0 nOutRF1]);
    
    
    sPlot_OutRF2Raster = subplot(5,2,9)
    imagesc(SmoothedSpikingMatrix(TrialIndexOutRF2(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesRaw nTimesRaw]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OutRF2Raster, 'Cue outside RF 2');
    ylabel('Trials');
    % ylim([0 nOutRF2]);
    xlabel('time (ms)');
    
    suptitle(['Neuron' int2str(n)]);
    set(f2, 'Position', [200 100 1000 1000]);
    
    print(f2,['Neuron ' int2str(n)], '-dpng');
    close;
end

%% plots of avg neuron firing rate and avg LFADS rate for in and off RF cue locations， plus spiking rasters
savedirTwo = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7'...
    '/TargetDim/postAnalysis/AvgFiringPlusRasterOnOffRFCueLoc/Gaussian50ms'];

cd(savedirTwo);

n = 1;
for n = 1:20
    f2 = figure
    
    sPlot_avgRealRates = subplot(3,2,1) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgFiringRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgFiringRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    
   
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 2
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgRealRates, 'Avg Real Firing Rate');
    ylabel('Firing Rate');
    set(gca,'xticklabel',{[]});
    
    xlabel('time (ms)');
%    legend('show')
    
    sPlot_avgLFADSRates = subplot(3,2,2) 
    % make the first subplot for plotting avg Reak rates for different cue locations
    plot(AvgLFADSRate_InRF(n,:), 'r', 'DisplayName','Cue in RF');
    hold on
    plot(AvgLFADSRate_OffRF(n,:), 'b', 'DisplayName','Cue opposite RF');
    
    
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_avgLFADSRates, 'Avg LFADS Firing Rate');
    ylabel('Firing Rate');
    set(gca,'YLim',sPlot_avgRealRates.YLim);
    xlabel('time (ms)');
    Legend = legend('show');
    set(Legend,'Position',[0.7387 0.5 0.1530 0.0655]);
    
    sPlot_InRFRaster = subplot(3,2,3)
    imagesc(WholeSpikingMatrix(TrialIndexInRF(n,:),:,n)) 
%     set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
%     set(gca,'XTickLabels',{'-800','TargetDim','+800'});
    set(gca,'xticklabel',{[]});
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_InRFRaster, 'Cue in RF');
    ylabel('Trials');
    % ylim([0 nInRF]);
    
    sPlot_OffRFRaster = subplot(3,2,5)
    imagesc(WholeSpikingMatrix(TrialIndexOffRF(n,:),:,n)) 
    set(gca,'XTick',[1 0.5*nTimesRaw nTimesRaw]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_OffRFRaster, 'Cue opposite RF');
    ylabel('Trials');
    % ylim([0 nOffRF]);
    xlabel('time (ms)');
    
    
    
    suptitle(['Neuron' int2str(n)]);
    set(f2, 'Position', [200 100 1000 900]);
    
    print(f2,['Neuron ' int2str(n)], '-dpng');
    close;
end
    
    


    
    

%% t-SNE analysis on rate vector over all trials, color and shape by cueLoc

% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_RelarrayOnsetCueLoc'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
alignType = 2; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType & holdIx == 0);
realAlign = r_realCopy.r(alignIx == alignType & holdIx == 0);

nTrials = length(lfadsAlign);
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
    title(['arrayOnset ' num2str(timeLag*4, '%+0.0f')]);
    print(f5,['arrayOnset ' num2str(timeLag*4, '%+0.0f')], '-dpng');
    close
end


%% t-SNE analysis on factor over all trials, color and shape by cueLoc

% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_onFactorArrayOnsetCueLoc'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
alignType = 2; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType & holdIx == 0);
realAlign = r_realCopy.r(alignIx == alignType & holdIx == 0);

nTrials = length(lfadsAlign);
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
    title(['arrayOnset ' num2str(timeLag*4, '%+0.0f')]);
    print(f5,['arrayOnset ' num2str(timeLag*4, '%+0.0f')], '-dpng');
    close
end


%%
%% t-SNE analysis on factor over all trials, color and shape by cueLoc

% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/' ...
    'CueOnArrayOnTargetDim_HoldRel/PostAnalysis/t-SNE_Analysis/t-SNE_onHeldoutSmoothedRateCueOnsetCueLoc'];

r_lfadsCopy = r_lfads(6).copy();
r_realCopy = r_real.copy();
alignType = 2; % select arrayOnset type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
lfadsAlign = r_lfadsCopy.r(alignIx == alignType & holdIx == 0);
realAlign = r_realCopy.r(alignIx == alignType & holdIx == 0);

nTrials = length(lfadsAlign);
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
    title(['arrayOnset ' num2str(timeLag*4, '%+0.0f')]);
    print(f5,['arrayOnset ' num2str(timeLag*4, '%+0.0f')], '-dpng');
    close
end


%%

% %% Colormap plot based on RT
% wcolormap = colormap('winter'); % 64 colors map in a sequence from blue to green
% conversionFactor = size(wcolormap, 1) / numel(rtVector);
% % convert the number of trials to the number of color, so every trial has a
% % color
% [~, trialOrder] = sort(rtVector, 'ascend'); 
% % sort trials based on RT. Each element in "trialOrder" indicates the index
% % of that element in original rtVector. In short, this means, if you sort
% % the rtVector based on RT, trialOrder tells you that, for this RT, which
% % trial it is.
% trialColors = zeros(numel(rtVector), 3);
% % initialize a trialColors matrix to tag each trial with a specific RT
% % color
% for ntrial = 1:numel(trialOrder)
%     colorIndex = ceil(ntrial * conversionFactor); 
%     % get the color index for that trial (the trial sorted in the sequence from fastest to slowest RT)
%     trialColors(trialOrder(ntrial), :) = wcolormap(colorIndex, :);
%     % use trialOrder to index to the original trial sequence and put the
%     % correct color there
% end
% 
% cd(savedirThree);
% f3 = figure
% % plot the tsne results and add coloring based on 'trialColors' matrix
% for ntrial = 1:nTrials
%     h = scatter(Rate_tsne(ntrial, 1), Rate_tsne(ntrial, 2));
%     set(h, 'markerfacecolor', trialColors(ntrial, :));
%     set(h, 'markeredgecolor','none');
%     set(h, 'Marker', trialShape(clVector(ntrial)));
%     % use the cued Location value in clVector to index to the proper trial
%     % shape in trialShape vector
%     hold on;
% end
% print(f3,['TargetDim+' int2str(timeLag)], '-dpng');
% close
% end


%% Heldout - predicted vs. real rt analysis

Predicted_rt = zeros(1, nTrials); 
% initiate a matrix to store predicted rt for each trial
RateMatrix_transform = RateMatrix';
% transform RateMatrix to prev Column x prev Row, this is good for
% alignment for calculate weights
TrialIndexRT = true(1, nTrials);
% Initialize a logical ones vector, to index which trial is held out for
% training

for nTrial = 1: nTrials
    
    TrialIndexRT(1,nTrial) = false;
    % put a logical zero to the trial that is held-out
    Trials_train = RateMatrix_transform(:, TrialIndexRT);
    % Get 302 trials for training
    rt_train = rtVector(TrialIndexRT);
    % get RTs for 302 trials for training
    w_NeuralToRT = mrdivide(rt_train, Trials_train);
    Predicted_rt(nTrial) = w_NeuralToRT * RateMatrix_transform(:,nTrial);
end
cd(savedirFive)
f5 = figure;
scatter(rtVector, Predicted_rt, 10, 'b','filled');
title(['TargetDim ' num2str(timeLag, '%+0.0f')]);
ylabel('Predicted RT');
xlabel('True RT');
print(f5,['TargetDim ' num2str(timeLag, '%+0.0f')], '-dpng');
close
%%

%% t-SNE analysis on rate vector over all trials, color and shape by cueLoc

