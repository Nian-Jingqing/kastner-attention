


%%
Day1_W = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withExternalInput_20180614/param_SGorjS/all/pbt_run/g061_w19/model_params'], '/LFADS_170127_cueOnArrayOnTargetDim_HoldRel.h5_out_fac_linear_W:0' );
%%
Day1_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withExternalInput_20180614/param_SGorjS/all/pbt_run/g061_w19/model_params'], '/LFADS_170127_cueOnArrayOnTargetDim_HoldRel.h5_out_fac_linear_b:0' );
%%

Day2_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170130_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );
%%
Day3_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170201_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day4_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170211_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day5_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170308_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day6_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170311_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day7_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170320_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day8_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170324_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day9_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170327_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day10_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170329_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day11_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170331_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

Day12_Wrates = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/' ...
                    'runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], '/LFADS_glm_fac_2_logrates_170407_cueOnArrayOnTargetDim_HoldRel.h5_W:0' );

%%
Day1_b = h5read(['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withGoodNeurons_Run_20180314/param_gnmwq8/all//lfadsOutput/model_params'], ...
                '/LFADS_glm_fac_2_logrates_170127_cueOnArrayOnTargetDim_HoldRel.h5_b:0' );
%%