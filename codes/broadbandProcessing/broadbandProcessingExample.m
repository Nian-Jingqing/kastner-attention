filename = '/mnt/scratch/feng/Remy_02262019_PUL_Raw.pl2';
outfile = '/mnt/scratch/cpandar/Remy_02262019_PUL_spikeband.mat';

tic;
broadband2streamMinMax( filename, outfile )
toc;

