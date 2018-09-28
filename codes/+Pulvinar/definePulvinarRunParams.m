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
par.c_batch_size = 210; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; 
par.useAlignmentMatrix = false; % use alignment matrices

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_ci_enc_dim = 64; % number of units in encoder for controller input
par.c_con_dim = 64; % number of units in controller
par.c_l2_gen_scale = 250;
par.c_l2_con_scale = 250;

par.c_learning_rate_stop = 1e-5; % we can stop really early for the demo

% Dataset parameters
par.spikeBinMs = 10; % data bin size
par.nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ]; % neuron indices

%##################### ADD PARAMETERS FOR RUN #######################
rc.addParams(par);
%####################################################################

%% RUN 2
%
% DESCRIPTION
%
% After the initial run using the combined dataset, we found some
% promising results using the factors and rates to decode the
% kinematic and EMG signals using single timepoint linear decoding.
% From comparing the linear decoding estimates from spiking data and
% factor data on opposing force field folds demonstrates that the
% factors seem to contain some information about the systme that
% allows for single timepoint decoding that the spiking data is
% unable to do. To determine the optimal number of factors to explain
% the underlying dynamics, I will sweep the factors hyperparameter to
% reduce the dimensionality of the factors and test its ability to
% perform single timepoint decoding on opposing force field folds of
% data. I will utilize this measure of prediciton fit as a means to
% determine the amount of information the factors are providing to
% perform optimal linear decoding. The sweep of the factors HP will
% consist of Runs 2-5.
%
%############ CONSERVED PARAMETERS FROM PREVIOUS RUN ################
% Model hyperparameters
par.c_co_dim = 3; % number of units in controller
par.c_batch_size = 210; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; 
par.useAlignmentMatrix = false; % use alignment matrices

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_ci_enc_dim = 64; % number of units in encoder for controller input
par.c_con_dim = 64; % number of units in controller
par.c_l2_gen_scale = 100;
par.c_l2_con_scale = 100;

par.c_learning_rate_stop = 1e-5; % we can stop really early for the demo

% Dataset parameters
par.spikeBinMs = 10; % data bin size

%############# ALTERED PARAMETERS FROM PREVIOUS RUN #################
par.nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ]; % neuron indices
%##################### ADD PARAMETERS FOR RUN #######################
rc.addParams(par);
%####################################################################







%% RUN 3
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
par.c_batch_size = 210; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; 
par.useAlignmentMatrix = false; % use alignment matrices

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_ci_enc_dim = 64; % number of units in encoder for controller input
par.c_con_dim = 64; % number of units in controller
par.c_l2_gen_scale = 250;
par.c_l2_con_scale = 250;

par.c_learning_rate_stop = 1e-6; % we can stop really early for the demo

% Dataset parameters
par.spikeBinMs = 10; % data bin size
par.nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ]; % neuron indices

%##################### ADD PARAMETERS FOR RUN #######################
rc.addParams(par);
%####################################################################

%% RUN 4
%
% DESCRIPTION
%
% After the initial run using the combined dataset, we found some
% promising results using the factors and rates to decode the
% kinematic and EMG signals using single timepoint linear decoding.
% From comparing the linear decoding estimates from spiking data and
% factor data on opposing force field folds demonstrates that the
% factors seem to contain some information about the systme that
% allows for single timepoint decoding that the spiking data is
% unable to do. To determine the optimal number of factors to explain
% the underlying dynamics, I will sweep the factors hyperparameter to
% reduce the dimensionality of the factors and test its ability to
% perform single timepoint decoding on opposing force field folds of
% data. I will utilize this measure of prediciton fit as a means to
% determine the amount of information the factors are providing to
% perform optimal linear decoding. The sweep of the factors HP will
% consist of Runs 2-5.
%
%############ CONSERVED PARAMETERS FROM PREVIOUS RUN ################
% Model hyperparameters
par.c_co_dim = 3; % number of units in controller
par.c_batch_size = 210; % must be < 1/5 of the min trial count
par.c_factors_dim = 20; 
par.useAlignmentMatrix = false; % use alignment matrices

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_ci_enc_dim = 64; % number of units in encoder for controller input
par.c_con_dim = 64; % number of units in controller
par.c_l2_gen_scale = 100;
par.c_l2_con_scale = 100;

par.c_learning_rate_stop = 1e-6; % we can stop really early for the demo

% Dataset parameters
par.spikeBinMs = 10; % data bin size

%############# ALTERED PARAMETERS FROM PREVIOUS RUN #################
par.nIndices = [ 8 10 13 32 38 40 41 51 52 60 70 71 81 86 96 102 ]; % neuron indices
%##################### ADD PARAMETERS FOR RUN #######################
rc.addParams(par);
%####################################################################










% %% RUN 3
% %
% % DESCRIPTION
% %
% % After the initial run using the combined dataset, we found some
% % promising results using the factors and rates to decode the
% % kinematic and EMG signals using single timepoint linear decoding.
% % From comparing the linear decoding estimates from spiking data and
% % factor data on opposing force field folds demonstrates that the
% % factors seem to contain some information about the systme that
% % allows for single timepoint decoding that the spiking data is
% % unable to do. To determine the optimal number of factors to explain
% % the underlying dynamics, I will sweep the factors hyperparameter to
% % reduce the dimensionality of the factors and test its ability to
% % perform single timepoint decoding on opposing force field folds of
% % data. I will utilize this measure of prediciton fit as a means to
% % determine the amount of information the factors are providing to
% % perform optimal linear decoding. The sweep of the factors HP will
% % consist of Runs 2-5.
% %
% %############ CONSERVED PARAMETERS FROM PREVIOUS RUN ################
% % Model hyperparameters
% par.c_co_dim = 3; % number of units in controller
% par.c_batch_size = 210; % must be < 1/5 of the min trial count
% par.c_factors_dim = 20; 
% par.useAlignmentMatrix = false; % use alignment matrices
% 
% par.c_gen_dim = 64; % number of units in generator RNN
% par.c_ic_enc_dim = 64; % number of units in encoder RNN
% par.c_ci_enc_dim = 64; % number of units in encoder for controller input
% par.c_con_dim = 64; % number of units in controller
% par.c_l2_gen_scale = 250;
% par.c_l2_con_scale = 250;
% 
% par.c_learning_rate_stop = 5e-5; % we can stop really early for the demo
% 
% % Dataset parameters
% par.spikeBinMs = 10; % data bin size
% 
% %############# ALTERED PARAMETERS FROM PREVIOUS RUN #################
% par.nIndices = [ 4 10 13 16 29 36 37 38 41 46 51 52 58 60 69 70 71 72 77 79 80 83 86 89 90 94 96 99 102 ]; % neuron indices
% %##################### ADD PARAMETERS FOR RUN #######################
% rc.addParams(par);
% %####################################################################
% 
% %% RUN 4
% %
% % DESCRIPTION
% %
% % After the initial run using the combined dataset, we found some
% % promising results using the factors and rates to decode the
% % kinematic and EMG signals using single timepoint linear decoding.
% % From comparing the linear decoding estimates from spiking data and
% % factor data on opposing force field folds demonstrates that the
% % factors seem to contain some information about the systme that
% % allows for single timepoint decoding that the spiking data is
% % unable to do. To determine the optimal number of factors to explain
% % the underlying dynamics, I will sweep the factors hyperparameter to
% % reduce the dimensionality of the factors and test its ability to
% % perform single timepoint decoding on opposing force field folds of
% % data. I will utilize this measure of prediciton fit as a means to
% % determine the amount of information the factors are providing to
% % perform optimal linear decoding. The sweep of the factors HP will
% % consist of Runs 2-5.
% %
% %############ CONSERVED PARAMETERS FROM PREVIOUS RUN ################
% % Model hyperparameters
% par.c_co_dim = 3; % number of units in controller
% par.c_batch_size = 210; % must be < 1/5 of the min trial count
% par.c_factors_dim = 20; 
% par.useAlignmentMatrix = false; % use alignment matrices
% 
% par.c_gen_dim = 64; % number of units in generator RNN
% par.c_ic_enc_dim = 64; % number of units in encoder RNN
% par.c_ci_enc_dim = 64; % number of units in encoder for controller input
% par.c_con_dim = 64; % number of units in controller
% par.c_l2_gen_scale = 100;
% par.c_l2_con_scale = 100;
% 
% par.c_learning_rate_stop = 5e-5; % we can stop really early for the demo
% 
% % Dataset parameters
% par.spikeBinMs = 10; % data bin size
% 
% %############# ALTERED PARAMETERS FROM PREVIOUS RUN #################
% par.nIndices = [ 4 10 13 16 29 36 37 38 41 46 51 52 58 60 69 70 71 72 77 79 80 83 86 89 90 94 96 99 102 ]; % neuron indices
% %##################### ADD PARAMETERS FOR RUN #######################
% rc.addParams(par);
% %####################################################################

%% ADD RUNSPECS

for iR = 1:dc.nDatasets
    runSpec = Pulvinar.RunSpec(dc.datasets(iR).getSingleRunName(), dc, dc.datasets(iR).name);
    rc.addRunSpec(runSpec);
end