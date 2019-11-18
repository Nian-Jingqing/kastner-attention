%% make pbt plots to compare the two runs: 180614 and 190508
addpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt/'))

%% plot pbt results for 190508
pbt_run_path = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190508/param_AsRGi7/all/pbt_run/';
results_dir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/reRun180614_20190508/param_AsRGi7/all/';
PBT_analysis.make_pbt_run_plots(pbt_run_path, results_dir)

%% plot pbt results for 180614
pbt_run_path_2 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withExternalInput_20180614/param_SGorjS/all/pbt_run/';
results_dir_2 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withExternalInput_20180614/param_SGorjS/all/';
PBT_analysis.make_pbt_run_plots(pbt_run_path_2, results_dir_2)