function rcs = defineMazeRunsetExample( runRoot )

%% path to maze data
dpathBase = '/snel/share/share/data/Shenoy/Maze/';

%% LOAD WELL CHARACTERIZED JENKINS 108 CONDITION DATASET
monkey = 'Jenkins';
dcj = Datasets.Maze.DatasetCollection( fullfile( dpathBase, monkey ) );
dcj.monkey = monkey;
dcj.autoDetectDatasets();

% the best-characterized Maze dataset is Jenkins 2009-09-18 [ dc.datasets(3) ]

%% LOAD THE NITSCHKE 108 CONDITION DATASETS
monkey = 'Nitschke';
dcn = Datasets.Maze.DatasetCollection( fullfile( dpathBase, monkey ) );
dcn.monkey = monkey;
dcn.autoDetectDatasets();

%  Nitschke datasets 2 & 3 both have 108 conditions



%%
% create an LFADS run collection for jenkins
rcj = Datasets.Maze.RunCollection(runRoot, 'mazeRunsJ', dcj);
rcj.version = 20180331;

%%
% create an LFADS run collection for nitschke
rcn = Datasets.Maze.RunCollection(runRoot, 'mazeRunsN', dcn);
rcn.version = 20180331;


%%
% set some hyperparams
par = Datasets.Maze.RunParams;
par.spikeBinMs = 10; % rebin the data at 2 ms
par.c_co_dim = 0; % no controller --> no inputs to generator
par.c_batch_size = 15; % must be < 1/5 of the min trial count

par.c_factors_dim = 10; % number of factors read out from generator to generate rates

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN

par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo

par.align = 'movementOnset';

%%
rcj.addParams( par );

%%
rcj.addRunSpec( Datasets.Maze.RunSpec( sprintf('d_%03i', 3 ), dcj, 3));

%%
rcn.addParams( par );

%%
rcn.addRunSpec( Datasets.Maze.RunSpec( sprintf('d_%03i', 3 ), dcn, 3));



rcs = { rcj rcn };