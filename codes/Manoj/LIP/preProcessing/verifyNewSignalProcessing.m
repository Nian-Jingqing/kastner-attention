%% addpath
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/CPspikepanels/utils');

%%
baseDir = '/mnt/scratch/feng/LIP';

dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
dataset(5).date = '03272019'; % previously didn't work
%dataset(5).date = '04062019';
dataset(6).date = '04252019';
dataset(7).date = '05022019';
outDir_new = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
outDir_old = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp1/';

%%
poolobj = gcp('nocreate');
delete( poolobj )
parpool( 20 )
for day = 5
    %day = 1; % let's verify the first session
    filename = ['Remy_RP_' dataset(day).date '_LIP_WB.pl2'];
    filedir = fullfile(baseDir, filename);
    outfilename = ['Remy_' dataset(day).date '_LIP_spikeband.mat'];
    outfiledir_new = fullfile(outDir_new, outfilename);
    outfiledir_old = fullfile(outDir_old, outfilename);

    %
    tic;
    broadband2streamMinMax( filedir, outfiledir_new, [300, 5000] );
    toc;
    clear functions
    
    %tic;
    %broadband2streamMinMax_oldSP( filedir, outfiledir_old );
    %toc;
end

%%

%%
%spikeBandFile = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/Remy_02182019_LIP_spikeband.mat';
for i = 7
    date = dataset(i).date;
    spikeBandBase = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
    spikeBandFileName = ['Remy_', date, '_LIP_spikeband.mat'];
    spikeBandFile = fullfile(spikeBandBase, spikeBandFileName);
    bb = load(spikeBandFile);
    bb = bb.spikeband;
    
    % plotting spiking panel
    plottingSpikePanel( bb, 'New Signal Processing Strategy', date)

    % --------------- start dealing with PSTH stuff ------------------%

    %
    plottingPSTHPanel(bb, date)
end

%%    
spikeBandFile = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/Remy_04062019_LIP_spikeband.mat';
bb = load(spikeBandFile);
bb = bb.spikeband;

%% plotting spiking panel
plottingSpikePanel( bb, 'New Signal Processing Strategy', '04062019')

% --------------- start dealing with PSTH stuff ------------------%

%%
%%
plottingPSTHPanel(bb, '04062019')
