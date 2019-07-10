% ----------- Preprocessing for running time warping ---------------
%% add paths
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/postAnalysisCodes')

%% set up datasets names
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

%% load alf
loadpath = '/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/timeWarping/';
dataName = 'spikes_rates_190620';
dataPath = fullfile(loadpath, dataName);
load(dataPath)
binsize_rescaled = 8;

%% load UEs
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

%%
% number of trials for each day
numTrialsTot = cellfun( @numel, alf );


% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial

minimalDelay = 700;
for nday = 1 : numel( alf )
    isCorrectArray{ nday } = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nday }.arrayShapesCorrect);
    isLongDelay{ nday } = ( [ alf{ nday }.arrayDim ] - [ alf{ nday }.arrayOnset] ) > ( minimalDelay/binsize_rescaled ); % try this later
    isCueLoc3{ nday } = UEs{ nday }.cueLoc == 3;
    trialsToKeep{ nday } = isCorrectArray{ nday } & ( isLongDelay{ nday } )' & isCueLoc3{ nday };% & UE2.isHoldTrial;
    %trialsToKeep{ nday } = isCueLoc3{ nday }; %for cueOnset jittering

    cueLocs{ nday } = unique(UEs{ nday }.cueLoc);

    for nc = 1 : numel( cueLocs{ nday } )
        trialsByCueLoc{ nday }{nc} = find( trialsToKeep{ nday } & (UEs{ nday }.cueLoc==cueLocs{ nday }(nc)));
        rtsByCueLoc{ nday }{nc} = UEs{ nday }.rt( trialsByCueLoc{ nday }{nc} );    
    end
end

%%%%%%%%%%%%%%%%%% from now, just do one day (day 5) for an example %%%%%%%%%%%%%%%%%%%%%%%

%% for cueOnset testing with jittering
% set up random jittering setting
nNeurons = size(alf{5}(1).spikes, 1);
trialsToKeepInds{ 5 } = find( trialsToKeep{ 5 } );
jitterBase = round(200*(rand(numel(trialsToKeepInds{ 5 }), 1) - 0.5)/binsize_rescaled); % create the amount of jittering to apply for each trial, uniform sampling (-100, 100)
baseWindow = round( [-100 300] / binsize_rescaled );
whichfieldTW = 'cueOnset';

% set up chopping window based on jittering on each trial and chop
newWindows = repmat(baseWindow, numel(jitterBase), 1) - jitterBase;
TW_tensor = zeros(numel(trialsToKeepInds{5}), numel(baseWindow(1):baseWindow(2)), nNeurons);
%TW_tensor = zeros(numel(trialsToKeepInds{5}), numel(baseWindow(1):baseWindow(2)), 1);
for itr = 1:numel(trialsToKeepInds{5})
    timePoints_tmp = newWindows(itr, 1):newWindows(itr, 2);
    ntr = trialsToKeepInds{ 5 }(itr);
    TW_tensor(itr, :, :) = alf{ 5 }( ntr ).spikes(:, alf{ 5 }( ntr ).( whichfieldTW ) + timePoints_tmp)';   
end

%% for arrayOnset
nNeurons = size(alf{5}(1).spikes, 1);
trialsToKeepInds{ 5 } = find( trialsToKeep{ 5 } );
baseWindow = round( [-100 300] / binsize_rescaled );
whichfieldTW = 'arrayOnset';

timePoints = baseWindow(1) : baseWindow(2);
TW_tensor = zeros(numel(trialsToKeepInds{5}), numel(timePoints), nNeurons);
for itr = 1:numel(trialsToKeepInds{5})
    ntr = trialsToKeepInds{ 5 }(itr);
    TW_tensor(itr, :, :) = alf{ 5 }( ntr ).spikes(:, alf{ 5 }( ntr ).( whichfieldTW ) + timePoints)';   
end


%% verify by plotting just one example neuron
example = squeeze(TW_tensor(:,:,19));
imagesc(example)

%%
dataDir = '/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/timeWarping/arrayOnset200_700';
[data_aligned, shifts] = TimeWarp.timeWarpAlign( dataDir, TW_tensor );

%% verify re-aligned data
example_2 = squeeze(data_aligned(:,:,19));
imagesc(example_2)
c = flipud(gray);
c_1 = [c(1:3,:); c(38:64,:)];
colormap(gca, c_1)

%%
%hist(double(shifts))

x = -120:0.1:120;
y = -120:0.1:120;

scatter(binsize_rescaled*jitterBase', binsize_rescaled*double(shifts));
hold on
plot(x,y)
xlabel('True shifts (ms)')
ylabel('Predicted shifts (ms)')
title('Day 170308')

%%
find(jitterBase == max(jitterBase))
