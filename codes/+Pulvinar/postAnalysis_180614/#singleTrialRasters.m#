%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools')
%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
    'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRelSepForAO_JanToMar';

% add individual datasets
Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170320_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170324_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170327_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170329_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170331_cueOnArrayOnTargetDim_HoldRel.mat');
% Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/runs'];
rc2 = Pulvinar.RunCollection(runRoot, 'withExternalInput_20180614', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc2.version = 20180614;

% this script defines the run params
% Pulvinar.multiDayDefinePulvinarRunParams;
% add the ones we want for this run
%for nrun = 1:numel( par4 )
%    rc2.addParams( par4( nrun ) );
%end

par = Pulvinar.RunParams;
par.spikeBinMs = 4;
par.c_batch_size = 2; % must be < 1/5 of the min trial count
par.trainToTestRatio = 50;
par.useAlignmentMatrix = true; % use alignment matrices

% the following params are specified in pbt script run manager
par.c_factors_dim = 30; 
par.c_in_factors_dim = 30;
% par.c_co_dim = 0; % number of units in controller
% par.c_gen_dim = 64; % number of units in generator RNN
% par.c_ic_enc_dim = 64; % number of units in encoder RNN
% par.c_ci_enc_dim = 64; % number of units in encoder for controller input
% par.c_con_dim = 64; % controller dimensionality
% par.c_l2_gen_scale = 10;
% par.c_l2_con_scale = 10;
% par.c_kl_ic_weight = 0.2;
% par.c_kl_co_weight = 0.2;
par.c_learning_rate_stop = 5e-4;
par.c_learning_rate_decay_factor = 0.95;
par.c_kl_increase_steps = 1;
par.c_l2_increase_steps = 1;
par.c_ext_input_dim = 6;
par.c_inject_ext_input_to_gen = true;
par.c_keep_ratio = 1;

% the following params make sure this is a pbt run
par.doPBT = true;
par.PBTscript = '/snel/home/fzhu23/Projects/Pulvinar/codes/+Pulvinar/drive_scripts/pbt_script_run_manager_180614.py';


rc2.addParams( par );
rc2.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));

%% change day and define Directories here
dayName = 170130;
dayStr = num2str(dayName);
day_id = 2;
r_id = 1;
lfpBasePath = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/continuousOverlapChop/';
realDataBasePath = '/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag'; % from the pre-processed spiking data
rawDataBasePath = '/snel/share/share/data/kastner/pulvinar/multi-unit/preAligned/data_raw/MarToJun/v10/';
lfpFilePattern = '_cueOnArrayOnTargetDim_HoldRel_lfp.mat';
realDataFilePattern = '_cueOnArrayOnTargetDim_HoldRel.mat';

%% load real data
realDataFileName = [dayStr, realDataFilePattern];
realDataFileDir = fullfile(realDataBasePath, realDataFileName);
olapChopped = load(realDataFileDir);
olapChopped = olapChopped.combinedData;

%% get real data
trueSpikes = [olapChopped.r.r.spikes];
nSpikes = size(trueSpikes, 2);
trial_time_ms = 500;
trial_olap_ms = 100;
out = olapChopped.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );

%% get LFADS data
run = rc2.runs(r_id);
r_lfads = olapChopped.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms);
assembled_lfads = R.Rstruct(r_lfads);

%%
tot_spikes = out.counts;
% nPieces = size(tot_spikes, 1);
%spikes = zeros(size(trueSpikes, 1), nSpikes);
spikes(length(olapChopped.r.r)).spikes = 1;
%% assembly
nChops = 0;
for r_trial = 1: length(olapChopped.r.r)
    spikes(r_trial).spikes = zeros(size(olapChopped.r.r(r_trial).spikes));
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

% %% trialize the continuous spikes
% r(length(olapChopped.r.r)).spikes = 1;
% past_nSpikes = 0;
% for i = 1: length(r)
%     nSpikesThisTrial = size(olapChopped.r.r(i).spikes,2);
%     start_ind = past_nSpikes +1;
%     end_ind = past_nSpikes + nSpikesThisTrial;
%     r(i).spikes = spikes(:, start_ind:end_ind);
%     past_nSpikes = end_ind;
% end
%%
assembled_real = R.Rstruct(spikes);

%% get timing info
subFolderName = ['M20', dayStr];
timingBaseDir = fullfile(rawDataBasePath, subFolderName, 'Gratings');
timingFilePattern = fullfile(timingBaseDir, '*.mat'); % Change to whatever pattern you need.
timingFiles = dir(timingFilePattern);
firstTimingFileName = timingFiles(1).name;
firstTimingFile = fullfile(timingBaseDir, firstTimingFileName);
data = load(firstTimingFile);
UE = data.UE;
clear data

%% get cueOn, arrayOn, targetDim timing from real data
rawSampleRate = 1000;
sessStartTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(1)*rawSampleRate;
% session start time - I set it to be the time when the monkey starts
% fixation in the first trial
sessEndTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(end)*rawSampleRate;
extraEndMs = 500;

startInds = round(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue * rawSampleRate - sessStartTime);
startInds(1) = 1;
stopInds = round(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice * rawSampleRate - sessStartTime);

cueInds = round(UE.cueOnset * rawSampleRate - sessStartTime);
arrayInds = round(UE.arrayOnset * rawSampleRate - sessStartTime);
dimInds = round(UE.targetDim * rawSampleRate - sessStartTime);

cueStart = cueInds - startInds;
arrayStart = arrayInds - startInds;
dimStart = dimInds - startInds(UE.isHoldTrial);
dimStart_real = dimStart;

stopTimes = stopInds(UE.isHoldTrial) - startInds(UE.isHoldTrial);

%% get binned timing for LFADS rate
cueStart_lfads = round(cueStart/par.spikeBinMs);
arrayStart_lfads = round(arrayStart/par.spikeBinMs);
dimStart_lfads = round(dimStart/par.spikeBinMs);


%% set up aligned trials for real data and lfads data

half_ms = 300;
olapChopped_hold = olapChopped.r.r(UE.isHoldTrial);
olapChopped_hold = olapChopped_hold((stopTimes - dimStart) > 302);

% for real data
assembled_hold_real = R.Rstruct(assembled_real.r(UE.isHoldTrial));
assembled_hold_real = R.Rstruct(assembled_hold_real.r((stopTimes - dimStart) > 302));
dimStart_real = dimStart_real((stopTimes - dimStart_real) > 302);

% for lfads data
assembled_hold_lfads = R.Rstruct(assembled_lfads.r(UE.isHoldTrial));
assembled_hold_lfads = R.Rstruct(assembled_hold_lfads.r((stopTimes - dimStart) > 302));

dimStart_lfads = dimStart_lfads((stopTimes - dimStart) > 302);


%% Align the trials for real data

% for cueAlign and arrayAlign
for i = 1:length(assembled_real.r)
    assembled_real.r(i).cueAlign = assembled_real.r(i).spikes(:, (cueStart(i)-(half_ms-1)):(cueStart(i)+half_ms) );
    assembled_real.r(i).arrayAlign = assembled_real.r(i).spikes(:, (arrayStart(i)-(half_ms-1)):(arrayStart(i)+half_ms) );
end

for i = 1:length(assembled_hold_real.r)
    assembled_hold_real.r(i).dimAlign = assembled_hold_real.r(i).spikes(:, (dimStart_real(i)-(half_ms-1)):(dimStart_real(i)+half_ms) );
end


%% align the trials for LFADS rate
half_lfads = 300/par.spikeBinMs;
for i = 1:length(assembled_lfads.r)
    assembled_lfads.r(i).cueAlign = assembled_lfads.r(i).rates(:, (cueStart_lfads(i)-(half_lfads-1)):(cueStart_lfads(i)+half_lfads) );
    assembled_lfads.r(i).arrayAlign = assembled_lfads.r(i).rates(:, (arrayStart_lfads(i)-(half_lfads-1)):(arrayStart_lfads(i)+half_lfads) );
end

for i = 1:length(assembled_hold_lfads.r)
    assembled_hold_lfads.r(i).dimAlign = assembled_hold_lfads.r(i).rates(:, (dimStart_lfads(i)-(half_lfads-1)):(dimStart_lfads(i)+half_lfads) );
end

%% get experiment info (nTrials, nTimes, nNeurons)

nTrials = length(assembled_real.r); % get trial number
nTimesRaw = size(assembled_real.r(1).cueAlign, 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(assembled_real.r(1).cueAlign, 1); % get neuron nubmer

%% assign condition type to the assembled struct
for i = 1: length(assembled_real.r)
    if UE.cueLoc(i) == 1
        assembled_real.r(i).cueCond = 1;
        if UE.isHoldTrial(i)
            assembled_real.r(i).arrayCond = 3;
        else
            assembled_real.r(i).arrayCond = 5;
        end
    elseif UE.cueLoc(i) == 3
        assembled_real.r(i).cueCond = 2;
        if UE.isHoldTrial(i)
            assembled_real.r(i).arrayCond = 4;
        else
            assembled_real.r(i).arrayCond = 6;
        end
    end
end

%% rebin real data
rbinned_cue = assembled_real.binData({'cueAlign'}, [par.spikeBinMs]);
rbinned_array = assembled_real.binData({'arrayAlign'}, [par.spikeBinMs]);


%% plotting the input with LFADS rate w/o real binned spiking
savedirRoot = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/singleTrialRasters/normalized/'];
savedirBase = fullfile(savedirRoot, dayStr);
condNames = {'CueOnsetCueLoc1', 'CueOnsetCueLoc3', 'ArrayOnsetCueLoc1Hold', 'ArrayOnsetCueLoc3Hold', 'ArrayOnsetCueLoc1Rel', 'ArrayOnsetCueLoc3Rel', 'TargetDimCueLoc1', 'TargetDimCueLoc3' };
% condIx = [r_lfadsWhole.r.conditionId];
% chopStart = 1*(nTimesLFADS/4) + 1;
% chopEnd = nTimesLFADS - 1*(nTimesLFADS/4);
% chopStart_lfp = 301;
% chopEnd_lfp = 1100;
% channelLabel = cell(1, size(R(1).lfp, 1));
% channelVector = size(R(1).lfp, 1) : -1 : 1;
% for c = 1:size(R(1).lfp, 1)
%     channelLabel{c} = int2str(channelVector(c));
% end
% 
% div = 0.5;
clear set

trialNum = 60;
nTimesLFADS = size(assembled_lfads.r(1).cueAlign, 2);
for t = 1: trialNum
    f1 = figure;
    
    % normalize the LFADS rates for each channel
    norm_lfads_cueAlign = normalize(assembled_lfads.r(t).cueAlign, 'centered');
    
    sp(1) = subplot(2, 1, 1);
    imagesc( norm_lfads_cueAlign );
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    set(gca,'XTickLabels',{'-300','cueOnset','+300'});
    title(sp(1), 'LFADS rates');
    ylabel('Multi-units');
%         set(sp(7), 'FontSize', 7);
%         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
%         sp(8) = subplot(7,2,14);
%         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);


    % plot raw spiking
    sp(2) = subplot(2, 1, 2);
    imagesc( rbinned_cue(t).cueAlign );
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    set(gca,'XTickLabels',{'-300','cueOnset','+300'});
    title(sp(2), 'Real Spiking');
    ylabel('Multi-units');
% %         set(sp(8), 'FontSize', 7);


%         set(sp(1), 'Position', [0.1300 0.5456 0.7750 0.6])
%         set(sp(2), 'Position', [0.1300 0.1539 0.7750 0.3119]);
    suptitle(condNames{assembled_real.r(t).cueCond});
%         set(f1, 'Position', [279 53 648 913]);
%         set(f1, 'Position', [375 67 1079 899]);
    savedirOne = fullfile(savedirBase, condNames{assembled_real.r(t).cueCond});
    cd(savedirOne)
    print(f1, ['Trial ' int2str(t)], '-dpng');
    close;
    
    
    f2 = figure;
    
    % normalize the LFADS rates for each channel
    norm_lfads_arrayAlign = normalize(assembled_lfads.r(t).arrayAlign, 'centered');
    
    sp_a(1) = subplot(2, 1, 1);
    imagesc( norm_lfads_arrayAlign );
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    set(gca,'XTickLabels',{'-300','arrayOnset','+300'});
    title(sp_a(1), 'LFADS rates');
    ylabel('Multi-units');
%         set(sp(7), 'FontSize', 7);
%         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
%         sp(8) = subplot(7,2,14);
%         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);


    % plot raw spiking
    sp_a(2) = subplot(2, 1, 2);
    imagesc( rbinned_array(t).arrayAlign );
    set(gca,'XTick',[1 0.5*nTimesLFADS nTimesLFADS]);
    set(gca,'XTickLabels',{'-300','arrayOnset','+300'});
    title(sp_a(2), 'Real Spiking');
    ylabel('Multi-units');
% %         set(sp(8), 'FontSize', 7);


%         set(sp(1), 'Position', [0.1300 0.5456 0.7750 0.6])
%         set(sp(2), 'Position', [0.1300 0.1539 0.7750 0.3119]);
    suptitle(condNames{assembled_real.r(t).arrayCond});
%         set(f1, 'Position', [279 53 648 913]);
%         set(f1, 'Position', [375 67 1079 899]);
    savedirTwo = fullfile(savedirBase, condNames{assembled_real.r(t).arrayCond});
    cd(savedirTwo)
    print(f2, ['Trial ' int2str(t)], '-dpng');
    close;
    
end




%%



for condType = 1:nCond
    trialsForThisCond = r_lfadsWhole.r(condIx == condType);
    trialsOfRealForThisCond = rbinned(condIx == condType);
    nTrialsThisCond = length(trialsForThisCond);
    R_selected = R(condIx == condType);
    if nTrialsThisCond > 40
        trialNum = 40;
    else
        trialNum = nTrialsThisCond;
    end
    for t = 1:trialNum
        lfadsRatesThisTrial = trialsForThisCond(t).rates;
        inputThisTrial = trialsForThisCond(t).controller_outputs;
        realSpikingThisTrial = trialsOfRealForThisCond(t).spikeCounts;
        f1 = figure;
%         for i = 1:nInputs
%             sp(i) = subplot(7, 2, i*2 - 1);
%             plot(inputThisTrial(i, chopStart:chopEnd), 'b');
%             set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
%             set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
%             set(sp(i), 'FontSize', 7);
%             xlim([0 div*nTimesLFADS]);
%             title(sp(i), ['Input ' int2str(i)]);
%         end
%         sp(7) = subplot(7,2,13);
%         sp(7) = subplot('position', [0.1300 0.0339 0.3347 0.150]);


        % plot lfp data
%         lfpThisTrial = R_selected(t).lfp_theta;
%         chopped_lfp = lfpThisTrial(:,chopStart_lfp:chopEnd_lfp)';
%         dis = 0.5*(max(chopped_lfp(:)) - min(chopped_lfp(:)));
%         base = size(chopped_lfp, 2):-1:1;
%         riseBy = base*dis;
%         scrolled_lfp = chopped_lfp + riseBy;
%         sp(1) = subplot(2, 1, 1);
%         plot(scrolled_lfp, 'b', 'LineWidth', 1);
%         set(gca,'YTick', fliplr(riseBy));
%         set(gca,'YTickLabels',channelLabel);
%         set(gca,'XTick',[1 div*0.5*nTimesLFADS*par7.spikeBinMs div*nTimesLFADS*par7.spikeBinMs]);
%         set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
%         title(sp(1), 'LFP (2 - 15 Hz)')
%         ylabel('channels')
        
        % plot LFADS rate
        sp(2) = subplot(2, 1, 1);
        imagesc( lfadsRatesThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(2), 'LFADS rates');
        ylabel('Multi-units');
%         set(sp(7), 'FontSize', 7);
%         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
%         sp(8) = subplot(7,2,14);
%         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);


        % plot raw spiking
        sp(3) = subplot(2, 1, 2);
        imagesc( realSpikingThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(3), 'Real Spiking');
        ylabel('Multi-units');
% %         set(sp(8), 'FontSize', 7);
        
        
%         set(sp(1), 'Position', [0.1300 0.5456 0.7750 0.6])
%         set(sp(2), 'Position', [0.1300 0.1539 0.7750 0.3119]);
        suptitle(condNames{condType});
%         set(f1, 'Position', [279 53 648 913]);
%         set(f1, 'Position', [375 67 1079 899]);
        savedirOne = fullfile(savedirBase, condNames{condType});
        cd(savedirOne)
        print(f1, ['Trial ' int2str(t)], '-dpng');
        close;
    end
end
        
  
%%





