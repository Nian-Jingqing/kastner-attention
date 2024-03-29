%%
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/bin/PBT_HP_opt/utils')
addpath('/snel/home/fzhu23/bin/PBT_HP_opt')

%%

pbt_dir = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withExternalInput_20180614/param_SGorjS/all/pbt_run/'];
results_save_dir = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/pbt_results/'];
image_format = '-dpng';

make_pbt_run_plots( pbt_dir, results_save_dir, image_format);