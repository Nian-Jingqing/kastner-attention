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
filename = ['Remy_RP_' dataset(i).date '_LIP_WB.pl2'];
filedir = fullfile(baseDir, filename);
outfilename = ['Remy_' dataset(i).date '_LIP_spikeband.mat'];
outfiledir = fullfile(outDir, outfilename);

tic;
bb = broadband2streamMinMax( filedir, outfiledir, [300, 5000] );
toc;

%% Remove NaN for minSpikeBand
for ich = 1:size(bb.minSpikeBand,2)
    whereNan_msb = find( isnan( bb.minSpikeBand( : , ich ) ) );
    chMsb{ ich } = bb.minSpikeBand(1:(whereNan_msb(1) - 1), ich);
end

%% Call CP's function to grab threshold crossings (NEED TO WORK ON THIS)

% should call the function calcThresholdCrossings. The function takes N x T data, but my data is in cell array. Need to check whether all channels have same length.
threshMultOrFixed = 4.5;
useMultiplier = true;
tStep = 1/40000;
windowLength = [];

[t,inds, wfstds, threshMultOrFixed] = calcThresholdCrossings(wf, threshMultOrFixed, windowLength, tStep, useMultiplier);

%% make plots (NEED TO WORK ON THIS)



