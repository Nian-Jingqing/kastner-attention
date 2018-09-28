%% ADJUST THESE DEFINITIONS FOR YOUR OWN PATHS TO:
% analysis_tools
addpath( '/snel/home/cpandar/analysis_tools/' );
% lfads-run-manager
addpath( '/snel/home/cpandar/lfads/lfads-run-manager/src/' );
% run_lfadslite.py
lpp = '/snel/home/cpandar/lfadslite/run_lfadslite.py';

%% setup the run collections
runRoot = '~/maze_runs/runs';
rcs = Examples.Maze.defineMazeRunsetExample( runRoot );
rcj = rcs{1}; % run for monkeyJ
rcn = rcs{2}; % run for monkeyN

% write out the data files
rcj.prepareForLFADS();
rcn.prepareForLFADS();

%% write the training shell scripts
lpp = '/snel/home/cpandar/lfadslite/run_lfadslite.py';
rcj.runs( 1 ).writeShellScriptLFADSTrain('cuda_visible_devices', 0, ...
                                         'display', 1, ...
                                         'lfadsPythonPath', lpp);
rcn.runs( 1 ).writeShellScriptLFADSTrain('cuda_visible_devices', 0, ...
                                         'display', 1, ...
                                         'lfadsPythonPath', lpp);

%% write the posterior mean shell scripts
rcj.runs( 1 ).writeShellScriptLFADSPosteriorMeanSample('cuda_visible_devices', 0, ...
                                                  'num_samples_posterior', 256, ...
                                                  'lfadsPythonPath', lpp);

