%% DEFINE PARAMETERS FOR ENTIRE RUN COLLECTION

par = Pulvinar.RunParams;

%####################################################################
%% RUN 1
%
% DESCRIPTION
%
% This is the intial run of the Cherian Combined Dataset using LFADS.
% The hope of using the combined dataset is to determine if providing
% LFADS with a rich unstructured dataset from two different force
% fields, we can force LFADS to find a lower dimensional rep. of the
% data by decoupling the effect of a single direction of the force
% field. By provding LFADS with data from both sessions, it cannot
% overfit the model to one set of data and must find underlying
% dynamics that explain both sets of data. We have chosen a 5 ms bin
% width for the data, using all the data from the dataset, 500 ms
% trial time w/ 100 ms overlap between trials.
%
%############## INITIAL PARAMETERS FOR LFADS RUN ####################
% Model hyperparameters
par.c_co_dim = 3; % number of units in controller
par.c_batch_size = 160; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; 
par.useAlignmentMatrix = true; % use alignment matrices

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_ci_enc_dim = 64; % number of units in encoder for controller input
par.c_con_dim = 64; % controller dimensionality
par.c_l2_gen_scale = 250;
par.c_l2_con_scale = 250;

par.c_learning_rate_stop = 1e-5; % we can stop really early for the demo

% Dataset parameters
par.spikeBinMs = 10; % data bin size
% par.nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ]; % neuron indices

%##################### ADD PARAMETERS FOR RUN #######################
% rc.addParams(par);
%####################################################################


%% CP & FZ, 2018-02-14
par2 = par.copy();

% decrease the step time for L2 and KL
par2.c_kl_increase_steps = 50;
par2.c_l2_increase_steps = 50;

% decrease the regularizers
par2.c_l2_gen_scale = 1;
par2.c_l2_con_scale = 1;
par2.c_kl_ic_weight = 0.1;
par2.c_kl_co_weight = 0.1;

% allow it a lot of factors
par2.c_factors_dim = 40; 

%rc.addParams(par2);

%% CP & FZ, 2018-02-16
set(1).l2 = [ 1 1 1 ];
set(1).kl = [0.2 0.5 0.8];
set(2).l2 = [ 10 50 100];
set(2).kl = [ 0.5 0.5 0.5];
runIndex = 0;
for nset = 1:2
    for nvalue = 1:3
        runIndex = runIndex + 1;
        parThisRun( runIndex ) = par2.copy();
        parThisRun( runIndex ).c_learning_rate_stop = 5e-4; % we can stop really early for the demo
        parThisRun( runIndex ).c_learning_rate_decay_factor = 0.95;
        parThisRun( runIndex ).c_l2_gen_scale = set( nset ).l2( nvalue );
        parThisRun( runIndex ).c_kl_ic_weight = set( nset ).kl( nvalue );
        parThisRun( runIndex ).c_kl_co_weight = set( nset ).kl( nvalue );
    end
end

%% FZ, 2018-02-18
par4( 1 ) = par2.copy();
par4( 1 ).c_learning_rate_stop = 5e-4;
par4( 1 ).c_learning_rate_decay_factor = 0.95;
par4( 1 ).c_l2_gen_scale = 1;
par4( 1 ).c_kl_ic_weight = 0.2;
par4( 1 ).c_kl_co_weight = 0.2;
par4( 1 ).c_co_dim = 6;

par4( 2 ) = par2.copy();
par4( 2 ).c_learning_rate_stop = 5e-4;
par4( 2 ).c_learning_rate_decay_factor = 0.95;
par4( 2 ).c_l2_gen_scale = 100;
par4( 2 ).c_kl_ic_weight = 0.5;
par4( 2 ).c_kl_co_weight = 0.5;
par4( 2 ).c_co_dim = 6;

%% FZ, 2018-03-12
par5( 1 ) = par2.copy();
par5( 1 ).c_learning_rate_stop = 5e-4;
par5( 1 ).c_learning_rate_decay_factor = 0.95;
par5( 1 ).c_l2_gen_scale = 1;
par5( 1 ).c_kl_ic_weight = 0.2;
par5( 1 ).c_kl_co_weight = 0.2;
par5( 1 ).c_co_dim = 4;%% ADD RUNSPECS

par5( 2 ) = par2.copy();
par5( 2 ).c_learning_rate_stop = 5e-4;
par5( 2 ).c_learning_rate_decay_factor = 0.95;
par5( 2 ).c_l2_gen_scale = 1;
par5( 2 ).c_kl_ic_weight = 0.2;
par5( 2 ).c_kl_co_weight = 0.2;
par5( 2 ).c_co_dim = 6;%% ADD RUNSPECS

par5( 3 ) = par2.copy();
par5( 3 ).c_learning_rate_stop = 5e-4;
par5( 3 ).c_learning_rate_decay_factor = 0.95;
par5( 3 ).c_l2_gen_scale = 1;
par5( 3 ).c_kl_ic_weight = 0.2;
par5( 3 ).c_kl_co_weight = 0.2;
par5( 3 ).c_co_dim = 8;%% ADD RUNSPECS

par5( 4 ) = par2.copy();
par5( 4 ).c_learning_rate_stop = 5e-4;
par5( 4 ).c_learning_rate_decay_factor = 0.95;
par5( 4 ).c_l2_gen_scale = 1;
par5( 4 ).c_kl_ic_weight = 0.2;
par5( 4 ).c_kl_co_weight = 0.2;
par5( 4 ).c_co_dim = 12;%% ADD RUNSPECS

%% FZ, 2018-03-14
par6( 1 ) = par5( 2 ).copy();
par6( 2 ) = par5( 4 ).copy();

%% FZ, 2018-03-29 PBT with factors_dim = 25
par7 = par6( 1 ).copy();
par7.c_factors_dim = 25;


% for iR = 1:dc.nDatasets
%     runSpec = Pulvinar.RunSpec(dc.datasets(iR).getSingleRunName(), dc, dc.datasets(iR).name);
%     rc.addRunSpec(runSpec);
% end

% rc.addRunSpec(Pulvinar.RunSpec('all', dc, 1:dc.nDatasets));