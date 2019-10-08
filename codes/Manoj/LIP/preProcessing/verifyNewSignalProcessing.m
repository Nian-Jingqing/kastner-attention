%% addpath
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/CPspikepanels/utils');

%%
baseDir = '/mnt/scratch/feng/LIP';
outDir = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
dataset(5).date = '03272019';
dataset(6).date = '04062019';
dataset(7).date = '04252019';
dataset(8).date = '05022019';

day = 1; % let's verify the first session
filename = ['Remy_RP_' dataset(day).date '_LIP_WB.pl2'];
filedir = fullfile(baseDir, filename);
outfilename = ['Remy_' dataset(day).date '_LIP_spikeband.mat'];
outfiledir = fullfile(outDir, outfilename);

%%
tic;
bb = broadband2streamMinMax( filedir, outfiledir, [300, 5000] );
toc;

%%
%spikeBandFile = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/Remy_02182019_LIP_spikeband.mat';
spikeBandFile = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/Remy_02182019_LIP_spikeband.mat';
bb = load(spikeBandFile);
bb = bb.spikeband;

%% Remove NaN for minSpikeBand
chVsb{32} = 0;
for ich = 1:size(bb.validSpikeBand,2)
    tic;
    tmp = bb.validSpikeBand(9062009 : 18124016, ich);
    %whereNan_vsb = find( isnan( bb.validSpikeBand( : , ich ) ) );
    whereNan_vsb = find( isnan( tmp ) )
    %chVsb{ ich } = bb.validSpikeBand(1:(whereNan_vsb(1) - 1), ich);
    if ~isempty(whereNan_vsb)
        chVsb{ ich } = tmp(1:(whereNan_vsb(1) - 1));
    else
        chVsb{ ich } = tmp;
    end
    toc;
end
%%
vsb_mat = [chVsb{:}];

%% Call CP's function to grab threshold crossings (NEED TO WORK ON THIS)

% should call the function calcThresholdCrossings. The function takes N x T data, but my data is in cell array. Need to check whether all channels have same length.
threshMultOrFixed = -4;
useMultiplier = true;
tStep = 1/40000;
windowLength = 40; % default is 30. Changed to 40 to match a 40000-equvalent sampling rate as default
wf = vsb_mat';

[t,inds, wfstds, threshMultOrFixed] = calcThresholdCrossings(wf, threshMultOrFixed, windowLength, tStep, useMultiplier);

%% make plots (NEED TO WORK ON THIS)
figure
for ich = 1:32
    subplot(4,8,ich)
    for iSpike = 1:numel(inds{ich})
        plot(vsb_mat((inds{ich}(iSpike)-10:inds{ich}(iSpike)+10), ich)')
        hold on
    end
    title(['Channel ' int2str(ich)])
    axis tight
end
suptitle('With new signal processing strategy')



%% do the same thing for old signal processing strategy
day = 1; % let's verify the first session
outDir = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp1/';
filename = ['Remy_RP_' dataset(day).date '_LIP_WB.pl2'];
filedir = fullfile(baseDir, filename);
outfilename = ['Remy_' dataset(day).date '_LIP_spikeband.mat'];
outfiledir = fullfile(outDir, outfilename);

%%
tic;
bb_1 = broadband2streamMinMax_oldSP( filedir, outfiledir );
toc;

%%
plottingSpikePanel( bb_1, 'Old Signal Processing Strategy')
