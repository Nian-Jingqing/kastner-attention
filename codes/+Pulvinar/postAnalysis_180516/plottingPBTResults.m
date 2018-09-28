%%
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/bin/PBT_HP_opt/utils')
addpath('/snel/home/fzhu23/bin/PBT_HP_opt')

%%

pbt_dir = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/lorenz_runs/runs/' ...
    'test_multiday_PBT_180517/PBT/param_mFB3v1/all/pbt_run/'];
results_save_dir = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'lorenz_runs/postAnalysis/test_multiday_PBT_180517/PBT/'];
image_format = '-dpng';

make_pbt_run_plots( pbt_dir, results_save_dir, image_format);