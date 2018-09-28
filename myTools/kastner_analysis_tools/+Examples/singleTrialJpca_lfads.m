%% set your paths
addpath('/snel/home/rbarker/Projects/matlab/analysis_tools/');

%% 
ddir = '/snel/share/share/derived/lfads_runs/Perich/perich_test2/run-output';
matname = 'lfads_run1_out.mat';
%matname = 'Chewie_20150319_CO_CS_BL_001_stripped.mat';

fname = fullfile( ddir, matname );
load(fname)
lfadsTimes = -130:2:249;
%% get some trial-averaged data
% movePrePost is written in terms of number of indices based on lfads window
movePrePost = [1 189]
avgLfads = run.alignAndAverage_lfads( 'conditionId', { 'rates', 'y_smoothed','spikes','kin' }, [1 189], lfadsTimes );


%% define some parameters for jPCA
softenNorm = 1;
suppressBWrosettes = true;
suppressHistograms = true;
numPCs = 10;
jPCA_params = struct();
jPCA_params.softenNorm = softenNorm;
jPCA_params.suppressBWrosettes = suppressBWrosettes;
jPCA_params.suppressHistograms = suppressHistograms;
jPCA_params.numPCs = numPCs;

%% Prepare for trial averaged jPCA
jPCA_times = -50:2:100;

% turn trial average data into input format for jpca 
avgJpcaFormatted = R.alignedStructToJpca( avgLfads, 'rates', 'conditionField', 'conditionCode') ;

%%  perform jpca on the trial-averaged data 
[ avgProjections, avgSummary ] = jPCA_with_CP_mods.jPCA( avgJpcaFormatted, jPCA_times, jPCA_params );

% put up some plots
avgColors = jPCA_with_CP_mods.phaseSpace( avgProjections, avgSummary );

%% Plot a few planes with avged data
plotParams.planes2plot = [1 2 3];
colors = jPCA_with_CP_mods.phaseSpace( avgProjections, avgSummary, plotParams)

%% plot extended times with same colors
plotParams.colors = avgColors.colors
for i = 1:numel(avgProjections)
    avgProjections(i).proj = avgProjections(i).projAllTimes;    
end
jPCA_with_CP_mods.phaseSpace( avgProjections,avgSummary, plotParams );

%% Prepare for single trial analysis
% convert this to a struct where each element is an individual trial 
%singleTrialPerCondition = R.singleTrialStructToAlignedStruct( run.r, 'rates' );
% convert this into input format for jpca
% We need a time vec, and we have previously defined it. For now, I will define it again here, but in the future, it should be defined in the lfads output code.
% The timeVecMs should be the same as any of the avged conditions
for i = 1:numel(run.r)
    run.r(i).timeVecMs = avgLfads(1).timeVecMs;
end

%% Single Trial jPCA
[ singleTrialProjections ] = jPCA_with_CP_mods.projectOntoJPCAPlane_lfads( avgJpcaFormatted, jPCA_times, jPCA_params, ...
                                                  avgSummary, run.r );
singleColors = jPCA_with_CP_mods.phaseSpace( singleTrialProjections, avgSummary );

%% Plot a few planes with Single Trial Data
plotParams.planes2plot = [1 2 3];
singleColors = jPCA_with_CP_mods.phaseSpace( singleTrialProjections, avgSummary, plotParams)

%% Plot Projection for all times with single trial data
%Specify colors to match the smaller time window
for i = 1:numel(colors)
    plotParams(i).colors = colors(i).colors;
end
%% Prepare to plot all times
for i = 1:numel(singleTrialProjections)
    singleTrialProjections(i).proj = singleTrialProjections(i).projAllTimes;
end

%% Make and save a movie with the data from all times (keeping colors the same)
movieParams.plane2plot = 1
movieParams.times = -130:2:249;
movieParams.colorsToUse =singleColors;
MV =jPCA_with_CP_mods.phaseMovie(singleTrialProjections, avgSummary, movieParams);
%% Save video to a file
% movie2avi has been removed from MATLAB
%movie2avi(MV, 'Perich_lfadsProjAllTimes_plane1.avi', 'FPS',12, 'compression','none')
% We will try using videowriter
mvObj = VideoWriter('Perich_lfadsProjAllTimes_plane1.avi');
mvObj.FrameRate = 12;
mvObj.open();
mvObj.writeVideo(MV);
mvObj.close();