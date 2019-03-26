%% add your paths here.

% add paths for Feng
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/postAnalysisCodes')

%% test jPCA code
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools/fromMarksLibraries')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools/CircStat2010d')

%% save the loaded data
%out.alf = alf;
%out.rfLoc = rfLoc;
%out.binsize_rescaled = binsize_rescaled;
%out.UEs = UEs;
saveDir2 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/loadedData/';
cd(saveDir2)
%saveName = 'allDays';
%save(saveName, 'out');
data = load('allDays');
alf = data.out.alf;
UEs = data.out.UEs;
binsize_rescaled = data.out.binsize_rescaled;

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

% also get rid of short delay trials
minimalDelay = 700;

for nday = 1 : numel( alf )
    isCorrectArray{ nday } = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nday }.arrayShapesCorrect);
    trialsToKeep{ nday } = isCorrectArray{ nday } & UEs{ nday }.cueLoc == 3;% & UE2.isHoldTrial;

    cueLocs{ nday } = unique(UEs{ nday }.cueLoc);

    for nc = 1 : numel( cueLocs{ nday } )
        trialsByCueLoc{ nday }{nc} = find( trialsToKeep{ nday } & (UEs{ nday }.cueLoc==cueLocs{ nday }(nc)));
        rtsByCueLoc{ nday }{nc} = UEs{ nday }.rt( trialsByCueLoc{ nday }{nc} );    
    end
end

%%

% concatenate all the factors

% % this window was used for finding oscillations during arrayDelay in
% %         factors 6 7 8 for session 6
window = round( [-300  0] / binsize_rescaled );



%window = round( [-1200 : 500] / binsize );
%window = round( [-200 1000] / binsize_rescaled );

 %window_raw = [-200 1200];

whichfieldDimred = 'arrayDim';
%whichfieldDimred = 'arrayOnset';

%window_lfp = round( [-296 32] / binsize_rescaled );

whichfieldJPCA = 'arrayDim';
window_jPCA = round( [-600  200] / binsize_rescaled );


%whichfieldPlot = 'arrayDim';
%newWindow = round( [-900  650] / binsize_rescaled );



whichfieldPlot = 'arrayOnset';
newWindow = round( [-200  1100] / binsize_rescaled );


%whichfieldPlot = 'cueOnset';
%newWindow = round( [-200  1200] / binsize_rescaled );


%% get factors into struct to prepare for jPCA
timePoints = window(1):window(2);
timePoints_jPCA = window_jPCA(1):window_jPCA(2);
%timePoints_lfp = window_lfp(1):window_lfp(2);
%timePoints_raw = window_raw(1):window_raw(2);
numBins = numel( timePoints);
numBins_jPCA = numel( timePoints_jPCA );
numFactors = size( alf{ nday }(1).rates, 1);
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );
%allFactors = zeros( numFactors, numBins * totalTrialsToKeep );

%%
for nday = 6
    ind = 1;
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        allFactors( :, (0:numBins-1) + ind ) = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        ind = ind + numBins;
    end
end

%%
%% do pca
meanFactors = mean( allFactors' );

[pca_proj_mat, pc_data] = pca( allFactors', 'NumComponents', 10);

clear dataForJPCA
dataForJPCA(totalTrialsToKeep).A = 0;
FRandSpiking = [];
ind = 1;
for nday = 1:numel( alf )
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        % dataForJPCA(ind).A = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        % dataForJPCA(ind).A = dataForJPCA(ind).A';
        frep = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        % mean center
        frep = frep - repmat( meanFactors(:), 1, numBins_jPCA );
        % project this data
        dim_reduced_data = pca_proj_mat' * frep;
        dataForJPCA(ind).A = dim_reduced_data(6:9,:)';
        dataForJPCA(ind).times = timePoints_jPCA*binsize_rescaled;
        dataForJPCA(ind).times = dataForJPCA(ind).times';
        ind = ind + 1;
    end
end

%% perform jPCA
jPCA_params.softenNorm = 5;
jPCA_params.suppressBWrosettes = true;
jPCA_params.suppressHistograms = true;
%%
times = -200:8:0;
jPCA_params.numPCs = 4;
[Projection, Summary] = jPCA(dataForJPCA, times, jPCA_params);
%%
plotParams.planes2plot = [ 1 2 ];
phaseSpace(Projection,Summary, plotParams);

%%
figure
trial_length = 12;
outdir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/traj_LowD/SFN/cueOnset/tmp8/';
cd(outdir);
cmap = lines();
for nt = 1 : 5 : numel(timePoints_jPCA)
%for nt = 71
    clf;
    for itr = 1:numel(Projection)
        startind = max(1, nt - trial_length);
        h = Plot.patchline(Projection(itr).projAllTimes(startind:nt, 1), Projection(itr).projAllTimes(startind:nt, 2), 'edgealpha',0.5);
        set( h, 'edgecolor', cmap(2,: ));
        set( h, 'facealpha', 0.5, 'edgealpha', 0.5 );
        hold on;

        h2 = plot(Projection(itr).projAllTimes(nt, 1), Projection(itr).projAllTimes(nt, 2), '.');
        set(h2, 'markersize', 5);
        set(h2, 'color', cmap(2,:) );
    end
    p = gca();
    axis(p, [ -0.15 0.1 -0.07 0.1]);
    title(timePoints_jPCA(nt)*binsize_rescaled);
    h_SFN = gcf();
    %    printpdf(h_SFN,int2str(nt) )
    pause(0.1);
    savefig(h_SFN, int2str(nt))
end
