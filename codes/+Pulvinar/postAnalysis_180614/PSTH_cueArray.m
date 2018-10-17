%%% add your paths here.

% add paths for Feng
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath(genpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes'))
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')

%%
buildRuns

 
%%
loadChoppedCombined
 

%% make a place to store output videos
outdir = '/Users/feng/SNEL/tmp/kastnervid/';
if ~isdir( outdir )
    mkdir( outdir );
end

%% only on laptop
% load data from file
%cp_paths_laptop
%load ~/tmp/forPlotting

%% fix any weirdness with zeros in the ALF
for nd = 1:6
    for ntr = 1:numel(alf{nd})
        alf{nd}(ntr).rates(alf{nd}(ntr).rates==0) = nan;
        alf{nd}(ntr).rt = UEs{nd}.rt( ntr );
    end
end


%%
% number of trials for each day
numTrialsTot = cellfun( @numel, alf );


% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial


for nday = 1 : numel( alf )
    isCorrectArray{ nday } = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nday }.arrayShapesCorrect);
    trialsToKeep{ nday } = isCorrectArray{ nday };% & UE2.isHoldTrial;
