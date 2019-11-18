%% addpath
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/CPspikepanels/utils');

%%

%dataset(1).date = '02182019';
%dataset(2).date = '03062019';
%dataset(3).date = '03112019';
%dataset(4).date = '03142019';
%dataset(5).date = '03272019'; % previously didn't work
%dataset(5).date = '04062019';
%dataset(6).date = '04252019';
%dataset(7).date = '05022019';


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

dataset(1).date = '02272019';
dataset(2).date = '03042019';
dataset(3).date = '03072019';
dataset(4).date = '03092019';
dataset(5).date = '03102019';
dataset(6).date = '03122019';
dataset(7).date = '03132019';
dataset(8).date = '03152019';

for i = 1:8
    date = dataset(i).date;
    spikeBandBase = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
    spikeBandFileName = ['Remy_', date, '_LIP_spikeband.mat'];
    spikeBandFile = fullfile(spikeBandBase, spikeBandFileName);
    bb = load(spikeBandFile);
    bb = bb.spikeband;
    for multiple = [3.5, 4, 4.5]        
    %for multiple = [6.5]        
        % plotting spiking panel
        tic;
        plottingSpikePanel( bb, 'New Signal Processing Strategy', date, multiple)

        % --------------- start dealing with PSTH stuff ------------------%

        %
        plottingPSTHPanel(bb, date, multiple)
        toc;
        clear functions
    end
end