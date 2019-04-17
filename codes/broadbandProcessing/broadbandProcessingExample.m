%filename = '/mnt/scratch/feng/Remy_02272019_PUL_RAW.pl2';
%outfile = '/mnt/scratch/feng/Remy_02272019_PUL_spikeband.mat';

%tic;
%broadband2streamMinMax( filename, outfile )
%toc;

%%
baseDir = '/mnt/scratch/feng';
dataset(1).date = '02082019';
dataset(2).date = '02182019';
dataset(3).date = '03082019';
dataset(4).date = '03102019';
dataset(5).date = '03112019';

for i = 1
    filename = ['Remy_RP_' dataset(i).date '_PUL_Raw.pl2'];
    filedir = fullfile(baseDir, filename);
    outfilename = ['Remy_' dataset(i).date '_PUL_spikeband.mat'];
    outfiledir = fullfile(baseDir, outfilename);

    tic;
    broadband2streamMinMax( filedir, outfiledir )
    toc;

end


    