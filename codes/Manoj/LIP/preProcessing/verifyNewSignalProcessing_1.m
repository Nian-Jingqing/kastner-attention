%% addpath
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/CPspikepanels/utils');

%%
w_id = 'finance:ftseries:tsmovavg:FunctionToBeRemoved';
warning('off',w_id)

%%
%baseDir = '/mnt/scratch/feng/LIP';
baseDir = '/home/feng';

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

%dataset(1).date = '02272019';
%dataset(2).date = '03042019';
%dataset(3).date = '03072019';
%dataset(4).date = '03092019';
%dataset(5).date = '03102019';
%dataset(6).date = '03122019';
%dataset(7).date = '03132019';
%dataset(8).date = '03152019';

dataset(1).date = '03162019';
dataset(2).date = '03182019';
dataset(3).date = '03292019';
dataset(4).date = '03312019';
dataset(5).date = '04012019';
dataset(6).date = '04032019';
dataset(7).date = '04052019';
dataset(8).date = '04242019';
dataset(9).date = '04262019';
dataset(10).date = '04292019';

outDir_new = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
outDir_old = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp1/';

%%
poolobj = gcp('nocreate');
delete( poolobj )
parpool( 12 )
for day = 1:numel(dataset)
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