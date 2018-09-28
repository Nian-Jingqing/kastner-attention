%% set your paths
addpath('/snel/home/rbarker/Projects/matlab/analysis_tools/');
addpath('/snel/home/rbarker/Projects/matlab/Perich/')
%% 
ddir = '/snel/share/share/data/Miller_Perich/';
matname = 'Chewie_20161011_CO_FF_BL_001_stripped.mat';
%matname = 'Chewie_20150319_CO_CS_BL_001_stripped.mat';

fname = fullfile( ddir, matname );

%%
[C, startInds, stopInds, trialstruct] = Perich_loadAndProcess( fname );

%%
% smooth the spiketrains with a gaussian kernel
C.smoothField( 'spikes', 'y_smoothed', 40 );

%%
clear r R
% turn into a trialized (R) struct
r = Movement.movementData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

%% trim to successful trials only
outcomes = [r.r.result];
r.r = r.r( outcomes == 'R' );


%% get some trial-averaged data
r.findMovementOnsetTimes(.1, 'rel');
movePrePost = [ 130 250 ];
avgGoAligned = r.alignAndAverage( 'conditionID', 'y_smoothed', movePrePost, 'moveOnset' );


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

%% Prepare for avg trial jPCA
jPCA_times = -50:100;
% turn trial average data into input format for jpca 
avgJpcaFormatted = R.alignedStructToJpca( avgGoAligned, 'y_smoothed', 'conditionField', 'conditionCode') ;
%  perform jpca on the trial-averaged data 
[ avgProjections, avgSummary ] = jPCA_with_CP_mods.jPCA( avgJpcaFormatted, jPCA_times, jPCA_params );

% put up some plots
avgColors = jPCA_with_CP_mods.phaseSpace( avgProjections, avgSummary );

%% Plot multiple planes from all times
% first establish colors from small time window
plotParams.planes2plot = [1 2 3];
colors = jPCA_with_CP_mods.phaseSpace(avgProjections, avgSummary, plotParams);

% Plot all Times
for i = 1:numel(colors)
    plotParams(i).colors = colors(i).colors;
end
% Define data for all times
for i = 1:numel(avgProjections)
    avgProjections(i).proj = avgProjections(i).projAllTimes;
end
jPCA_with_CP_mods.phaseSpace(avgProjections, avgSummary, plotParams)

%% get the same data as above, but individual trials
singleTrialGoAligned = r.getAligned( 'conditionID', 'y_smoothed', movePrePost, 'moveOnset' );
% convert this to a struct where each element is an individual trial 
singleTrialPerCondition = R.singleTrialStructToAlignedStruct( singleTrialGoAligned, 'y_smoothed' );
% convert this into input format for jpca 
singleTrialJpcaFormatted = R.alignedStructToJpca( singleTrialPerCondition, 'y_smoothed', 'conditionField', 'conditionCode') ;

%% Run jPCA on single trial data
[ singleTrialProjections ] = jPCA_with_CP_mods.projectOntoJPCAPlane( avgJpcaFormatted, jPCA_times, jPCA_params, ...
                                                  avgSummary, singleTrialJpcaFormatted );
jPCA_with_CP_mods.phaseSpace( singleTrialProjections, avgSummary );

%% Prepare to make single trial movie
singleColors = jPCA_with_CP_mods.phaseSpace( singleTrialProjections, avgSummary)
for i = 1:numel(singleTrialProjections)
    singleTrialProjections(i).proj = singleTrialProjections(i).projAllTimes;
end


%% Make and save a movie with the data from all times (keeping colors the same)
movieParams.plane2plot = 1
movieParams.times = -130:249;
% We only need to plot one plane
movieParams.colorsToUse =singleColors(1);
MV =jPCA_with_CP_mods.phaseMovie(singleTrialProjections, avgSummary, movieParams);

%% Save video to a file
% movie2avi has been removed from MATLAB
%movie2avi(MV, 'Perich_lfadsProjAllTimes_plane1.avi', 'FPS',12, 'compression','none')
% We will try using videowriter
mvObj = VideoWriter('Perich_lfadsProjAllTimes_plane1.avi');
mvObj.FrameRate = 24;
mvObj.open();
mvObj.writeVideo(MV);
mvObj.close();