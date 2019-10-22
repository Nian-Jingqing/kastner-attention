%% addpath
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/CPspikepanels/utils');

%%
%baseDir = '/mnt/scratch/feng/LIP';

dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
%dataset(5).date = '03272019';
dataset(5).date = '04062019';
dataset(6).date = '04252019';
dataset(7).date = '05022019';
thresh_multiple = [4, 4, 4, 4, 4, 4, 5];

%%
for i = 3:7
    date = dataset(i).date;
    spikeBandBase = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
    spikeBandFileName = ['Remy_', date, '_LIP_spikeband.mat'];
    spikeBandFile = fullfile(spikeBandBase, spikeBandFileName);
    bb = load(spikeBandFile);
    bb = bb.spikeband;
    spikeband_to_spikes_func(bb, date, thresh_multiple(i))
    clear functions
end