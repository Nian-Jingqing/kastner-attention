%% plot costs - non-pbt
addpath('/snel/home/fzhu23/bin/pbt-hp-opt')
%addpath('/snel/home/mreza/projects/lfads_utils')

f1 = figure;
%fitlog_name = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/runs/tc_day0218_noNan_190822/param_kJIY08/all/lfadsOutput/fitlog.csv';
%fitlog_name = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/runs/tc_7Days_noNan_190822/param_BsX8Qs/all/lfadsOutput/fitlog.csv';
fitlog_name = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/runs/tc_7Days_withNan_190822/param_BsX8Qs/all/lfadsOutput/fitlog.csv';
[val,col] = PBT_analysis.read_fitlog(fitlog_name);
rec_train = val(:,col(['recon' '_train']));
rec_train = rec_train(~isnan(rec_train));
plot(rec_train, 'DisplayName', ['lfadslite' ' ' 'recon' ' train'])
hold on
rec_valid = val(:,col(['recon' '_valid']));
rec_valid = rec_valid(~isnan(rec_valid));
plot(rec_valid, 'DisplayName', ['lfadslite' ' ' 'recon' ' valid'])
ylabel('Loss')
legend('-DynamicLegend');
min_loss = rec_valid(end);
title(['M1 run (final recon loss valid = ' num2str(min_loss) ')'])
set(f1, 'Position', [190 53 1411 886])



%% plot costs - pbt
%% make pbt plots
%rmpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt/'));
%addpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt_modifiedPlotting_190607/'));
%rmpath(genpath('/snel/home/fzhu23/bin/PBT_HP_opt_modifiedPlotting_190607/'));
addpath(genpath('/snel/home/fzhu23/bin/pbt-hp-opt/'));

%% plot pbt results for the run
%pbt_run_path = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/runs/tc_PBT_190822/param_8k7EL7/all/pbt_run/';
pbt_run_path = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/runs/tc_noNan_PBT_10msBin_190826/param_NGvjDx/all/pbt_run/';
%pbt_run_path = '/snel/share/share/derived/kastner/model_interaction/runs/PMd_190809/param_f8zU1j/all/pbt_run/';
results_dir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/PBT_10msBin_noDataMask_20190826/costs/';
%results_dir = '/snel/share/share/derived/kastner/model_interaction/postAnalysis/maze_analysis/PMd_run190809/pbt_plots/';
PBT_analysis.make_pbt_run_plots(pbt_run_path, results_dir)