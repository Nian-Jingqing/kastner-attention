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

%%
load exampleData
jPCA_params.softenNorm = 5;
jPCA_params.suppressBWrosettes = true;
jPCA_params.suppressHistograms = true;
%%
times = 50:10:300;
jPCA_params.numPCs = 6;
[Projection, Summary] = jPCA(Data, times, jPCA_params);
%%
phaseSpace(Projection,Summary);
printFigs(gcf, '.', '-dpdf', 'Basic jPCA plot');

%%
buildRuns_20180614

%%
loadChoppedCombined_twoLocations
 
%% load and preprocess LFP
loadLFP_twoLocations

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
    isLongDelay{ nday } = ( [ alf{ nday }.arrayOnset ] - [ alf{ nday }.cueOnset] ) > ( minimalDelay/binsize_rescaled );
    trialsToKeep{ nday } = isCorrectArray{ nday } & (~[lfp{ nday }.isTrialOutlier])' & ( isLongDelay{ nday } )';

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

window_lfp = round( [400 700] / binsize_rescaled );

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
timePoints_lfp = window_lfp(1):window_lfp(2);
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


%%
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

%%
clear FRandSpiking
FRandSpiking = [];
ind = 1;
for nday = 1:numel( alf )
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        FRandSpiking(ind).FR = alf{ nday }( ntr ).FR( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        FRandSpiking(ind).spiking = alf{ nday }( ntr ).spikes( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        ind = ind + 1;
    end
end

%% pick up correct trials, filter, resample and chop LFP
filtHighCutoff = 8;
filtLowCutoff = 2;
Fs = 1000;
chopResampled_lfp(totalTrialsToKeep).lfps = 0;
ind = 1;
for nday = 1:numel( lfp )
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        tmp_filteredLFP = bandpassFilter_singleTrial( lfp{ nday }( ntr ).lfps, filtHighCutoff, filtLowCutoff, Fs);
        tmp_resampledLFP = (resample(tmp_filteredLFP', 1,binsize_rescaled))';
        chopResampled_lfp(ind).lfps = tmp_resampledLFP( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints_lfp );
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
f1 = figure;
for i = 1:10
    subplot(5,2,i);
    plot(Projection(i+499).proj(:,1), 'r');
    hold on
    plot(Projection(i+499).proj(:,2), 'b');
    hold on
end


%%
savedir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/jPCA/';
cd(savedir)
printFigs(gcf, '.', '-dpdf', '400To500ms_plane1');

%% get projections from jPCA
%dimId = 1;
%proj_rates(numel(Projection)). = 0;
%for nTrial = 1 : numel(Projection)
    %proj_rates(nTrial) = Projection(nTrial).proj(:, dimId)';
    %end


%% plotting single trials plots


%%
trialsToKeepEachDay = cellfun( @sum, trialsToKeep );
%%
maxLag = 80;
ind = 0;
%crossCorr{ numel( alf ) } = [];
crossCorr = {};
for nday = 1:numel( alf )
    trialIds = ( ind + 1 ):( ind + trialsToKeepEachDay( nday ));
    for n = 1 : size(lfp{ nday }(1).lfps, 1)
        for nTrial = 1 : numel( trialIds )
            for dimId = 1 : 4
                fieldToStore = ['dim' num2str(dimId)];
                crossCorr{ nday }( n ).( fieldToStore )( nTrial, :) = xcorr(Projection( trialIds( nTrial ) ).proj( :, dimId)', chopResampled_lfp( trialIds( nTrial ) ).lfps( n, :), maxLag/binsize_rescaled);
            end
        end
    end
    ind = ind + trialsToKeepEachDay( nday );
end

                
%% plotting cross-correlation

savedir2_base = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/jPCA_LFP/tmp/';

clear set
timeLagStr = num2str(maxLag);
for nday = 1 : numel( crossCorr )
    savedir2 = [savedir2_base datasets( nday ).shortName];
    for n = 1 : numel( crossCorr{ nday } )
        f1 = figure;
        sp1 = subplot(2,2,1);
        shadedErrorBar([], crossCorr{ nday }( n ).dim1, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim1, 1)) }, 'lineProps', '-r')
        set(gca,'XTick',[1 0.5*size(crossCorr{ nday }( n ).dim1, 2) size(crossCorr{ nday }( n ).dim1, 2)]);
        timeStrLeft = ['-', timeLagStr];
        timeStrRight = ['+', timeLagStr];
        set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
        set(gca,'XLim',[0 size(crossCorr{ nday }( n ).dim1, 2)]);
        ylabel('cross-correlation');
        title('Dim 1');
        hold on
        
        sp2 = subplot(2,2,2);
        shadedErrorBar([], crossCorr{ nday }( n ).dim2, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim2, 1)) }, 'lineProps', '-r')
        set(gca,'XTick',[1 0.5*size(crossCorr{ nday }( n ).dim2, 2) size(crossCorr{ nday }( n ).dim2, 2)]);
        timeStrLeft = ['-', timeLagStr];
        timeStrRight = ['+', timeLagStr];
        set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
        set(gca,'XLim',[0 size(crossCorr{ nday }( n ).dim2, 2)]);
        ylabel('cross-correlation');
        title('Dim 2');
        hold on

        sp3 = subplot(2,2,3)
        shadedErrorBar([], crossCorr{ nday }( n ).dim3, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim3, 1)) }, 'lineProps', '-r')
        set(gca,'XTick',[1 0.5*size(crossCorr{ nday }( n ).dim3, 2) size(crossCorr{ nday }( n ).dim3, 2)]);
        timeStrLeft = ['-', timeLagStr];
        timeStrRight = ['+', timeLagStr];
        set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
        set(gca,'XLim',[0 size(crossCorr{ nday }( n ).dim3, 2)]);
        ylabel('cross-correlation');
        title('Dim 3');
        hold on

        sp4 = subplot(2,2,4)
        shadedErrorBar([], crossCorr{ nday }( n ).dim4, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim4, 1)) }, 'lineProps', '-r')
        set(gca,'XTick',[1 0.5*size(crossCorr{ nday }( n ).dim4, 2) size(crossCorr{ nday }( n ).dim4, 2)]);
        timeStrLeft = ['-', timeLagStr];
        timeStrRight = ['+', timeLagStr];
        set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
        set(gca,'XLim',[0 size(crossCorr{ nday }( n ).dim4, 2)]);
        ylabel('cross-correlation');
        title('Dim 4');
        
        suptitle(['Cross Correlation for LFP channel ' int2str(n)]);
        set(f1, 'Position', [229 79 1573 887]);
        cd(savedir2)
        print(f1,['Channel ' int2str(n)], '-dpng');
        close
    end
end