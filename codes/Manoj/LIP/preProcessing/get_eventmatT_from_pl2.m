addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/Manoj/LIP/preProcessing/')

%% set pl2 file path
baseDir = '/mnt/scratch/feng/LIP';
outDir = '/snel/share/share/data/kastner/Manoj/LIP/eventMat/';
dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
dataset(5).date = '03272019';
dataset(6).date = '04062019';
dataset(7).date = '04252019';
dataset(8).date = '05022019';

%%
for i = 1:8
    filename = ['Remy_RP_' dataset(i).date '_LIP_WB.pl2'];
    filedir = fullfile(baseDir, filename);
    outfilename = ['eventmat_' dataset(i).date '.mat'];
    %outfiledir = fullfile(baseDir, outfilename);

    tic;
    readEventCodes( filedir, outDir, outfilename );
    toc;

end