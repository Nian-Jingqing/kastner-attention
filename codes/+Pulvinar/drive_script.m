%% build the dataset collection
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
datasetPath = '/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/TargetDim/datasets';



%% Locate and specify the datasets
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'TargetDimOneSess';

% add individual datasets
Pulvinar.Dataset(dc, 'TargetDim_holdTrials.mat');
% MyExperiment.Dataset(dc, 'dataset002.mat');%Example code for lorenz
% example - for future use
% MyExperiment.Dataset(dc, 'dataset003.mat');

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = '/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/TargetDim/runs';
rc = Pulvinar.RunCollection(runRoot, 'exampleRun', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc.version = 20171212;

%% Set parameters for the entire run collection

par = Pulvinar.RunParams;
par.spikeBinMs = 2; % rebin the data at 2 ms
par.c_co_dim = 0; % no controller --> no inputs to generator
par.c_batch_size = 60; % must be < 1/5 of the min trial count % total trial number is 303, so I chose this to be 60
par.c_factors_dim = 10; % and manually set it for multisession stitched models
par.useAlignmentMatrix = false; % use alignment matrices initial guess for multisession stitching

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN

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

%% Prepare LFADS input

% generate all of the data files LFADS needs to run everything
rc.prepareForLFADS();

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
%rc.writeShellScriptRunQueue('display', 50, 'maxTasksSimultaneously', 4, 'gpuList', [0 1], 'virtualenv', 'tensorflow');
rc.writeShellScriptRunQueue('display', 9, 'maxTasksSimultaneously', 4, 'gpuList', [0 1]);


%% Post-running analysis - loading data and the output of LFADS
r_id = 1;
r = rc.runs(r_id); % pull out run information
r.loadSequenceData(); % load sequence data in that run
r.loadPosteriorMeans(); % load posterior mean in that run
r.addPosteriorMeansToSeq();
rate_struct = r.sequenceData{1}; % Put sequence data into a struct
nTrials = size(rate_struct,2); % get trial number
nTimesRaw = size(rate_struct(1).y, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(rate_struct(1).y, 1); % get neuron nubmer
nTimesLFADS = size(rate_struct(1).rates,2);% get trial length for rebinned data that was operated by LFADS model
data = dc.datasets(1).loadData(); % get the original dataset
%% Post-running analysis - Plots of LFADS rates and spiking activity for example trials

savedirOne = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7/TargetDim/' ...
'postAnalysis/LFADSRatesVsSpikingActivity'];

cd(savedirOne);
for i = 1:size(rate_struct,2) % loop over all trials
    f1 = figure
    sPlot_LFADSRates = subplot(2,1,1) % make the first subplot for plotting LFADS rates (all neurons) for the trial
    imagesc(rate_struct(i).rates) % "rate_struct(i).y" contains the LFADS predicted rates of all neurons for the trial
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_LFADSRates, 'LFADS rates');
    ylabel('neurons')
%     xlabel('time (ms)'); % delete this x lable to make the figure look
%     not that compact
    hold on
    plot([0.5*nTimesLFADS 0.5*nTimesLFADS],[1 nNeurons],'g'); % Plot a vertical line at aligned time
    
    
    sPlot_Spiking = subplot(2,1,2) % make the second subplot for plotting spiking activity (all neurons) for the trial
    imagesc(rate_struct(i).y); % "rate_struct(i).y" contains the spiking activity of all neurons for the trial
    set(gca,'XTick',[1 0.5*nTimesRaw nTimesRaw]); 
    % find the three points that need to be labeled on the x axis (start time, TargetDim time, and end time)
    set(gca,'XTickLabels',{'-800','TargetDim','+800'}); 
    % Put the appropriate labels on the x axis. * remember to change this
    % if necessary
    title(sPlot_Spiking, 'Spiking Activity');
    ylabel('neurons');
    xlabel('time (ms)');
    hold on
    plot([0.5*nTimesRaw 0.5*nTimesRaw],[1 nNeurons],'g'); % Plot a vertical line at aligned time
    
    suptitle(['Trial' int2str(i)]);
    
    print(f1,['Trial ' int2str(i)], '-dpng');
    close;
end

%% Post-running analysis - Compare avg neuron firing rate and avg LFADS rate for different cue locations
% prepare for plotting


WholeSpikingMatrix = permute(cat(3, rate_struct.y), [3 2 1]);
% re-organize the spiking data to nTrials x nTimes x nNeurons
WholeLFADSRatesMatrix = permute(cat(3, rate_struct.rates), [3 2 1]);
% re-organize the LFADS rates to nTrials x nTimes x nNeurons


clVector = arrayfun(@(x) x.cueLoc, data.R); 
% pull out cue location for each trial and put into a vector
rfLoc = data.R(1).rfloc;
% pull out receptive field for each neuron
AllCueLoc = unique(clVector);
rebinSize = 2; % data was rebinned from 1ms to 2ms during LFADS
SmoothedSpikingMatrix = zeros(nTrials, nTimesLFADS, nNeurons);
% initialize a spiking matrix for all neurons, with size nTrials x
% n rebinnedTimes x nNeurons
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

for n = 1:nNeuron % loop over all neurons

ReshapeStore = reshape(WholeSpikingMatrix(:,:,n), [nTrials rebinSize nTimesLFADS]);
% reshape the spiking data for neuron n to be able to rebin the data
RebinnedSpikingData = squeeze(sum(ReshapeStore, 2));
% Sum up spike counts in each bin to get the rebinned spiking data for this
% neuron

SmoothedSpikingMatrix(:,:,n) = (gaussianSmooth(RebinnedSpikingData', 50, rebinSize, 0))';
% apply Gaussian smoothing to the spiking data for each neuron. Gaussian
% window is 10 ms here. Change to 20ms or 50 ms later
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
savedirTwo = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7' ...
    '/TargetDim/postAnalysis/AvgFiringPlusRasterCueLocs/Gaussian10ms'];

cd(savedirTwo);

n = 1;
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
    
    sPlot_OffRFRaster = subplot(5,2,5)
    imagesc(WholeSpikingMatrix(TrialIndexOffRF(n,:),:,n)) 
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
    imagesc(WholeSpikingMatrix(TrialIndexOutRF1(n,:),:,n)) 
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
    imagesc(WholeSpikingMatrix(TrialIndexOutRF2(n,:),:,n)) 
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
    
    


    
    

%% t-SNE analysis on rate vector over all trials, color by reaction time
savedirThree = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7' ...
    '/TargetDim/postAnalysis/t-sne_Analysis/tsne_BasedOnColorMap'];
% path to store t-sne plot based on color map

savedirFour = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7' ...
    '/TargetDim/postAnalysis/t-sne_Analysis/FS_separateByMedian+-0.05'];
% path to store t-sne plot based on coarse separation
savedirFive = ['/snel/home/fzhu23/Projects/Pulvinar/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v7' ...
    '/TargetDim/postAnalysis/Held-Out_rtPrediction'];


timeLagAll = -25:5:25;
for i = 1:length(timeLagAll)


    timeLag = timeLagAll(i);
    % set a time lag relative to TargetDim for selecting neuron data at a
    % specific time point for analysis;
    TimeTargetDim = 400; 
    % the time at TargetDim is 800 ms, so it's the 400th in LFADS rates
    TimeForAnalysis = TimeTargetDim + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    rtVector = arrayfun(@(x) x.rt, data.R);
    % put all the reaction time into a vector
    clVector = arrayfun(@(x) x.cueLoc, data.R); 
    % pull out cue location for each trial and put into a vector
    
    for i = 1:nTrials
        RateMatrix(i,:) = rate_struct(i).rates(:,TimeForAnalysis);
        % pull out the neural data for all neuron at a specific time in that trial,
        % and put them into the ith row in the rate matrix
    end
    Rate_tsne = tsne(RateMatrix);
    % Actually run tsne
    trialShape = ['o', '>', 's', 'd'];
    % Use different shapes to indicate cue location of the trial.
    % circle for cueLoc = 1; right-pointing triangle for cueLoc = 2;
    % square for cueLoc = 3; diamond for cueLoc = 4;


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
%% Color plot based on coarse seperation based on RT
% 
%  rtTagSlow = rtVector > median(rtVector)+0.05;
%   % any rt > median rt + 0.05s in rtVector is considered to be slow, indicated by
%  % logical 1
%  rtTagFast = rtVector < median(rtVector)-0.05;
%  % any rt < median rt - 0.05s in rtVector is considered to be fast, indicated by
%  % logical 1
% 
% 
% rtSlowRates = RateMatrix(rtTagSlow,:);
% % Use rtTag to indice trials with slow rt
% rtFastRates = RateMatrix(rtTagFast,:);
% % Use rtTag to indice trials with fast rt
% Rate_tsne = tsne(RateMatrix);
% % apply tsne to reduce the dimensionality of 106 neurons to 2d
% cd(savedirFour)
% f4 = figure
% scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
% hold on
% scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
% title(['TargetDim+' int2str(timeLag)]);
% ylabel('2nd Dimension');
% xlabel('1st Dimension');
% legend('show')
% print(f4,['TargetDim+' int2str(timeLag)], '-dpng');
% close
end

%%




r_id = 1;
r = rc.runs(r_id); % pull out run information
r.loadSequenceData(); % load sequence data in that run
r.loadPosteriorMeans(); % load posterior mean in that run
r.addPosteriorMeansToSeq();
rate_struct = r.sequenceData{1}; % Put sequence data into a struct
nTrials = size(rate_struct,2); % get trial number
f = [rate_struct.factors]; % take out the factor matrix
f = [f; ones(1,size(f,2))]; % add the bias to factor matrix
data = dc.datasets(1).loadData(); % load the dataset
s_original = cat(2, data.R.spikeCounts); % re-organize the dataset to nChannels x n(Time*trial)
new_sR = size(f,2); % set new sample rate for the system
orig_sR = size(s_original,2); % get old sample rate from the original system
s = resample(s_original',new_sR,orig_sR)'; % resample original dataset to get to same length of f
train_Size = 0.7*size(s,2); % get the training group size for applying affine transformation
s_Train = s(:,1:train_Size); % get the training group of the data
f_Train = f(:,1:train_Size); % get the training group of the factors
s_Test = s(:,train_Size+1:end); % get the testing group of the data
f_Test = f(:,train_Size+1:end); % get the testing group of the factors
W_trans = mrdivide(s_Train,f_Train); % perform affine transformation
s_Predict = W_trans*f_Test; % get the predicted system
get_r_value = corrcoef(s_Test,s_Predict); 
r_square = get_r_value(1,2)^2; % calculate r square value to see how well the factors can predict the system
r_square

% %% Try something else
% 
% r_id = 1;
% r = rc.runs(r_id); % pull out run information
% r.loadSequenceData(); % load sequence data in that run
% r.loadPosteriorMeans(); % load posterior mean in that run
% r.addPosteriorMeansToSeq();
% rate_struct = r.sequenceData{1}; % Put sequence data into a struct
% nTrials = size(rate_struct,2); % get trial number
% f = [rate_struct.factors]; % take out the factor matrix
% f = [f; ones(1,size(f,2))]; % add the bias to factor matrix
% f = f(:, 1:121200);
% data = dc.datasets(1).loadData(); % load the dataset
% s_original = cat(2, data.R.spikeCounts); % re-organize the dataset to nChannels x n(Time*trial)
% new_sR = 242400; % set new sample rate for the system
% orig_sR = size(s_original,2); % get old sample rate from the original system
% s = resample(s_original',new_sR,orig_sR)'; % resample original dataset to get to same length of f
% s = s(:,1:121200);
% train_Size = 0.7*size(s,2); % get the training group size for applying affine transformation
% s_Train = s(:,1:train_Size); % get the training group of the data
% f_Train = f(:,1:train_Size); % get the training group of the factors
% s_Test = s(:,train_Size+1:end); % get the testing group of the data
% f_Test = f(:,train_Size+1:end); % get the testing group of the factors
% W_trans = mrdivide(s_Train,f_Train); % perform affine transformation
% s_Predict = W_trans*f_Test; % get the predicted system
% get_r_value = corrcoef(s_Test,s_Predict); 
% r_square = get_r_value(1,2)^2; % calculate r square value to see how well the factors can predict the system
%  r_square


%% Post-running analysis for lorenz system 

% r_id = 1;
% r = rc.runs(r_id);
% r.loadSequenceData();
% r.loadPosteriorMeans();
% r.addPosteriorMeansToSeq();
% 
% rate_struct = r.sequenceData{1};
% 
% nTrials = size(rate_struct,2);
% condID = [rate_struct.conditionId];
% nCond = length(unique(condID));
% nTpC = nTrials/nCond;
% 
% f = [rate_struct.factors];
% f = [f; ones(1,size(f,2))];
% 
% data1 = dc.datasets(1).loadData();
% 
% l_traj = data1.lorenz_trajectories;
% new_sR = 500;
% orig_sR = 1000;
% s = zeros(size(l_traj,1),size(f,2));
% total_count1 = 1;
% total_count2 = nTpC*new_sR;
% for i = 1:nCond
%     tmp = l_traj(:,:,i);
%     tmp_resamp = resample(tmp',new_sR,orig_sR)';
%     s(:,total_count1:total_count2) = repmat(tmp_resamp,1,nTpC);
%     total_count1 = total_count2+1;
%     total_count2 = total_count2+nTpC*new_sR;
% end
% W_trans = mrdivide(s,f);
%     
% %%
% 
% r_id2 = 2;
% r2 = rc.runs(r_id);
% r2.loadSequenceData();
% r2.loadPosteriorMeans();
% r2.addPosteriorMeansToSeq();
% 
% rate_struct2 = r2.sequenceData{1};
% 
% 
% 
% f2 = [rate_struct2.factors];
% f2 = [f2; ones(1,size(f,2))];
% 
% data2 = dc.datasets(2).loadData();
% 
% l_traj2 = data2.lorenz_trajectories;
% s2_actual = zeros(size(l_traj2,1),size(f2,2));
% total_count1 = 1;
% total_count2 = nTpC*new_sR;
% for i = 1:nCond
%     tmp2 = l_traj2(:,:,i);
%     tmp_resamp2 = resample(tmp2',new_sR,orig_sR)';
%     s2_actual(:,total_count1:total_count2) = repmat(tmp_resamp2,1,nTpC);
%     total_count1 = total_count2+1;
%     total_count2 = total_count2+nTpC*new_sR;
% end
% 
% s2_predicted = W_trans*f2;
% %%
% 
% r = corrcoef(s2_actual,s2_predicted);
% r_square = r(1,2)^2;
% r_square