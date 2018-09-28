classdef RunParams < LFADS.RunParams
   properties
       align = 'MoveOnsetOnline';
       nTrialsKeep = 0;
       
       % for maze data
       whichChannels = []
       whichTrials = []
       startOffset = 0;
       endOffset = 0;

       c_do_train_encoder_only = false;
       c_do_reset_learning_rate = false;
       c_checkpoint_pb_load_name = '';

       paramSetName = ''
   end
   
   methods
       function p = RunParams()
           % set default values here
           p.spikeBinMs = 2;
           p.whichChannels = [];
           p.startOffset = -450;
           p.endOffset = 450;

           p.whichTrials = [];

           % command line args
           p.c_batch_size = 128; % trials
           p.c_learning_rate_decay_factor = 0.95;
           p.c_l2_increase_steps = 2000;
           p.c_kl_increase_steps = 2000;
           p.c_keep_prob = 0.95; % dropout of units in the network
           p.c_gen_dim = 100; % generator network size
           p.c_ci_enc_dim = 128; % network size for controller input encoder
           p.c_ic_enc_dim = 128; % network size for IC encoder
           p.c_con_dim = 128; %controller dimensionality
           p.c_factors_dim = 30;
           p.c_co_dim = 4;
           p.c_ic_dim = 64;
           p.c_in_factors_dim = 50; % i don't beieve this is used in non-multiday 

           p.c_allow_gpu_growth = true;
           p.c_do_train_encoder_only = false;
           p.c_do_reset_learning_rate = false;
           p.c_checkpoint_pb_load_name = '';

       end
       
       function suffix = generateSuffix(p)
           if isempty(p.paramSetName)
               error('you must name your parameter set');
           end
           suffix = p.paramSetName;
       end
   end
end