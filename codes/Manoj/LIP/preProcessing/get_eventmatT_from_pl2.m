addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/Manoj/LIP/preProcessing/')

%% set pl2 file path
baseDir = '/mnt/scratch/feng/LIP';
outDir = '/snel/share/share/data/kastner/Manoj/LIP/eventMat/';

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

%dataset(1).date = '02182019';
%dataset(2).date = '03062019';
%dataset(3).date = '03112019';
%dataset(4).date = '03142019';
%dataset(5).date = '03272019';
%dataset(6).date = '04062019';
%dataset(7).date = '04252019';
%dataset(8).date = '05022019';

%dataset(1).date = '02272019';
%dataset(2).date = '03042019';
%dataset(3).date = '03072019';
%dataset(4).date = '03092019';
%dataset(5).date = '03102019';
%dataset(6).date = '03122019';
%dataset(7).date = '03132019';
%dataset(8).date = '03152019';

%dataset(1).date = '03162019';
%dataset(2).date = '03182019';
%dataset(3).date = '03292019';
%dataset(4).date = '03312019';
%dataset(5).date = '04012019';
%dataset(6).date = '04032019';
%dataset(7).date = '04052019';
%dataset(8).date = '04242019';
%dataset(9).date = '04262019';
%dataset(10).date = '04292019';

dataset(1).date = '03082019';

%%
for i = numel(dataset)
    filename = ['Remy_RP_' dataset(i).date '_LIP_WB.pl2'];
    filedir = fullfile(baseDir, filename);
    outfilename = ['eventmat_' dataset(i).date '.mat'];
    %outfiledir = fullfile(baseDir, outfilename);

    tic;
    readEventCodes( filedir, outDir, outfilename );
    toc;

end