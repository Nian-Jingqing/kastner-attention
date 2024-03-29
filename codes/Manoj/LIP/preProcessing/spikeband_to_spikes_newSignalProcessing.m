%% addpath
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/CPspikepanels/utils');

%%
%baseDir = '/mnt/scratch/feng/LIP';

% below are selected sessions for LFADS run, some sessions with lower threshold
%dataset(1).date = '02082019';
%dataset(2).date = '02132019';
%dataset(3).date = '02142019';
%dataset(4).date = '02152019';
%dataset(5).date = '02182019';
%dataset(6).date = '02262019';
%dataset(7).date = '02282019';
%dataset(8).date = '03032019';
%dataset(9).date = '03062019';
%dataset(10).date = '03142019';

% below are the new sessions that Manoj sent Oct 30th
%dataset(1).date = '02082019';
%dataset(2).date = '02132019';
%dataset(3).date = '02142019';
%dataset(4).date = '02152019';
%dataset(5).date = '02162019';
%dataset(6).date = '02262019';
%dataset(7).date = '02282019';
%dataset(8).date = '03012019';
%dataset(9).date = '03022019';
%dataset(10).date = '03032019';

% below are the first sessions that Manoj sent
%dataset(1).date = '02182019';
%dataset(2).date = '03062019';
%dataset(3).date = '03112019';
%dataset(4).date = '03142019';
%dataset(5).date = '03272019';
%dataset(5).date = '04062019';
%dataset(6).date = '04252019';
%dataset(7).date = '05022019';

% below are selected sessions from the 3rd and 4th cohort of sessions that Manoj sent
dataset(1).date = '03162019';
dataset(2).date = '03312019';
dataset(3).date = '04012019';
dataset(4).date = '04052019';
dataset(5).date = '04262019';

%thresh_multiple = [4, 4, 4, 4, 4, 4, 5]; % thresholds for first 7 sessions Manoj sent
%thresh_multiple = [3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 4.5]; % lower the threshold
%thresh_multiple = [4, 4, 4, 4, 4.5, 6, 4, 4, 4, 4]; % for new sessions (02082019 - 03032019)
%thresh_multiple = [3.5, 3.5, 3.5, 4, 3.5, 6, 4, 4, 3.5, 3.5]; % for lowering thresholds for selected sessions for LFADS run % works the best so far before 11/18/2019 (for 1st + 2nd cohort)
thresh_multiple = [3.5, 3.5, 3.5, 3.5, 3.5]; % for 3rd and 4th cohort

%%
for i = 1:numel(dataset)
    tic;
    date = dataset(i).date;
    %spikeBandBase = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
    spikeBandBase = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/4th_cohort_noLowPass/';
    spikeBandFileName = ['Remy_', date, '_LIP_spikeband.mat'];
    spikeBandFile = fullfile(spikeBandBase, spikeBandFileName);
    bb = load(spikeBandFile);
    bb = bb.spikeband;
    spikeband_to_spikes_func(bb, date, thresh_multiple(i))
    clear functions
    toc;
end