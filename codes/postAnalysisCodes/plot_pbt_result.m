%% make pbt plots to compare the two runs: 180614 and 190508
rmpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt/'));
addpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt_modifiedPlotting_190607/'));
%rmpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt_modifiedPlotting_190607/')); addpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt/'));

%% plot pbt results for 190508
pbt_run_path = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190508/param_AsRGi7/all/pbt_run/';
results_dir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190508/param_AsRGi7/all/';
PBT_analysis.make_pbt_run_plots(pbt_run_path, results_dir)

%% plot pbt results for 180530
pbt_run_path_2 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190530/param_mfUfYK/all/pbt_run/';
results_dir_2 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190530/param_mfUfYK/all/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_2, results_dir_2)

%% plot pbt results for 190603
pbt_run_path_2 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190603/param_U9qrDy/all/pbt_run/';
results_dir_2 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190603/param_U9qrDy/all/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_2, results_dir_2)

%% plot pbt results for 190605
pbt_run_path_3 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190605/param_WVw2ln/all/pbt_run/';
results_dir_3 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190605/param_WVw2ln/all/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_3, results_dir_3)

%% plot pbt results for 190607
pbt_run_path_4 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190607/param_BRr_3C/all/pbt_run/';
results_dir_4 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190607/param_BRr_3C/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_4, results_dir_4)

%% plot pbt results for 190615
pbt_run_path_5 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190615/param_LA14-f/all/pbt_run/';
results_dir_5 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190615/param_LA14-f/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_5, results_dir_5)

%% plot pbt results for 190616
pbt_run_path_6 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190616/param_Tk11H5/all/pbt_run/';
results_dir_6 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190616/param_Tk11H5/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_6, results_dir_6)

%% plot pbt results for 190618
pbt_run_path_7 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190618/param_U5c32w/all/pbt_run/';
results_dir_7 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190618/param_U5c32w/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_7, results_dir_7)

%% plot pbt results for 190619
pbt_run_path_8 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190619/param_gip6h5/all/pbt_run/';
results_dir_8 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190619/param_gip6h5/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_8, results_dir_8)

%% plot pbt results for 190620
pbt_run_path_9 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190620/param_YEGdZA/all/pbt_run/';
results_dir_9 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190620/param_YEGdZA/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_9, results_dir_9)

%% plot pbt results for 190622
pbt_run_path_10 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190622/param_yvX31y/all/pbt_run/';
results_dir_10 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190622/param_yvX31y/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_10, results_dir_10)

%% plot pbt results for 190623
pbt_run_path_11 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190623/param_Z8wEMU/all/pbt_run/';
results_dir_11 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190623/param_Z8wEMU/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_11, results_dir_11)

%% plot pbt results for 190626
pbt_run_path_12 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190626/param__sTJ-Z/all/pbt_run/';
results_dir_12 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190626/param__sTJ-Z/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_12, results_dir_12)

%% plot pbt results for 190627
pbt_run_path_13 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190627/param_pSXbYQ/all/pbt_run/';
results_dir_13 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190627/param_pSXbYQ/all/tmp_smoothed/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_13, results_dir_13)

%% plot the PBT result for rossler 
pbt_run_path_14 = '/snel/share/share/derived/kastner/model_interaction/runs/interactingRossler_200timesteps_191006/downstreamRossler_191006/pbt_run/'
results_dir_14 = '/snel/share/share/derived/kastner/model_interaction/runs/interactingRossler_200timesteps_191006/downstreamRossler_191006/'
PBT_analysis.make_pbt_run_plots(pbt_run_path_14, results_dir_14)

%% plot the PBT result for rossler 
pbt_run_path_15 = '/snel/share/share/derived/kastner/model_interaction/runs/downstreamLorenz_191018/pbt_run/'
results_dir_15 = '/snel/share/share/derived/kastner/model_interaction/runs/downstreamLorenz_191018/hp_progression/'
PBT_analysis.make_pbt_run_plots(pbt_run_path_15, results_dir_15)