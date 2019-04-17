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


%% load and preprocess LFP
%loadLFP_twoLocations



outdir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/pulvinar/180208/postAnalysis/withExternalInput_20190315/lowD_traj/jPCA/tmp/';
if ~isdir(outdir)
    mkdir(outdir);
end

cd(outdir)



%%
% number of trials for each day
numTrialsTot = cellfun( @numel, alf );


% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial

for nday = 1 : numel( alf )
    trialsToKeep{ nday } = (UE{ nday }.barType == 1 & ismember(UE{ nday }.cueType, [4]));
    trialsToKeep{ nday } = trialsToKeep{nday} & (([alf{nday}.targetStart] - [alf{nday}.cueOnset])*binsize_rescaled > 800);
    errorsToKeep{ nday } = trialsToKeep{ nday } & UE{ nday }.isErrorTrial & ~UE{ nday }.isEarlyError;
    hitsToKeep{ nday } = trialsToKeep{ nday } & ~UE{ nday }.isErrorTrial;
    cueLocs{ nday } = unique(UE{ nday }.cueType);
    

    for nc = 1 : numel( cueLocs{ nday } )
        trialsByCueLoc{ nday }{nc} = find( trialsToKeep{ nday } & (UE{ nday }.cueType==cueLocs{ nday }(nc)));
        %        rtsByCueLoc{ nday }{nc} = UEs{ nday }.rt( trialsByCueLoc{ nday }{nc} );    
    end
end

%%

% concatenate all the factors

% % this window was used for finding oscillations during arrayDelay in
% %         factors 6 7 8 for session 6
window = round( [400  800] / binsize_rescaled );



%window = round( [-1200 : 500] / binsize );
%window = round( [-200 1000] / binsize_rescaled );

 %window_raw = [-200 1200];

whichfieldDimred = 'cueOnset';
%whichfieldDimred = 'arrayOnset';

%window_lfp = round( [-296 32] / binsize_rescaled );

whichfieldJPCA = 'cueOnset';
window_jPCA = round( [0  800] / binsize_rescaled );


%whichfieldPlot = 'arrayDim';
%newWindow = round( [-900  650] / binsize_rescaled );



whichfieldPlot = 'cueOnset';
newWindow = round( [0  800] / binsize_rescaled );


%whichfieldPlot = 'cueOnset';
%newWindow = round( [-200  1200] / binsize_rescaled );


%% dimred based on all days
timePoints = window(1):window(2);
timePoints_jPCA = window_jPCA(1):window_jPCA(2);
numBins = numel( timePoints );
numBins_jPCA = numel( timePoints_jPCA );
numFactors = size( alf{ nday }(1).factors, 1);
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );
allFactors = zeros( numFactors, numBins * totalTrialsToKeep );

%
ind = 1;
for nday = 1 : numel( alf)
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        allFactors( :, (0:numBins-1) + ind ) = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        ind = ind + numBins;
    end
end 


%% calculate number of neurons for each day
for nday = 1 : numel( alf )
    nNeurons(nday) = size(alf{nday}(1).rates, 1);
end

%% load all weight matrix mapping from factors to LFADS rates
weights = [];
rootPath = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/pulvinar/180208/runs/withExternalInput_20190315/param_DghsNv/all/pbt_run/g058_w28/model_params';
for nday =1: numel(alf)
    dayFile_W = ['/LFADS_' datasets(nday).shortName '_v1.h5_out_fac_linear_W:0'];
    dayFile_b = ['/LFADS_' datasets(nday).shortName '_v1.h5_out_fac_linear_b:0'];
    weights(nday).W = h5read([rootPath], dayFile_W);
    weights(nday).W_T = weights(nday).W';
    weights(nday).b = h5read([rootPath], dayFile_b);
    weights(nday).b_T = weights(nday).b';
end

%% concatenate all dahys' weights together
Wrates_W_cat = [weights.W_T]';
Wrates_b_cat = [weights.b_T]';

%% project from factor space to rate space
allRates = Wrates_W_cat*(allFactors) + Wrates_b_cat;



%% do pca
meanRates = mean( allRates' );

[pca_proj_mat, pc_data] = pca( allRates', 'NumComponents', 30);

%% do pca - factor space
meanFactors = mean( allFactors' );

[pca_proj_mat, pc_data] = pca( allFactors', 'NumComponents', 10);


%%
clear dataForJPCA
dataForJPCA(totalTrialsToKeep).A = 0;
FRandSpiking = [];
ind = 1;
for nday = 1:numel( alf )
%for nday = [1 3]
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        % dataForJPCA(ind).A = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        % dataForJPCA(ind).A = dataForJPCA(ind).A';
        frep = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        %frep = Wrates_W_cat*(frep) + Wrates_b_cat; % if use rate space
        % mean center
        %frep = frep - repmat( meanRates(:), 1, numBins_jPCA ); % if use rate space
        frep = frep - repmat( meanFactors(:), 1, numBins_jPCA ); % if use factor space
        % project this data
        dim_reduced_data = pca_proj_mat' * frep;
        dataForJPCA(ind).A = dim_reduced_data(1:10,:)';
        dataForJPCA(ind).times = timePoints_jPCA*binsize_rescaled;
        dataForJPCA(ind).times = dataForJPCA(ind).times';
        ind = ind + 1;
    end
end

%% if do jPCA without doing PCA at all
clear dataForJPCA
dataForJPCA(totalTrialsToKeep).A = 0;
FRandSpiking = [];
ind = 1;
%RT = zeros(1, totalTrialsToKeep);
%DL = zeros(1, totalTrialsToKeep);
for nday = 1:numel( alf )
%for nday = [1 3]
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        dataForJPCA(ind).A = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        dataForJPCA(ind).A = dataForJPCA(ind).A';        
        dataForJPCA(ind).times = timePoints_jPCA*binsize_rescaled;
        dataForJPCA(ind).times = dataForJPCA(ind).times';
        %        RT(ind) = alf{ nday }( ntr ).rt;
        %        DL(ind) = alf{ nday }( ntr ).arrayDim - alf{ nday }( ntr ).arrayOnset;
        ind = ind + 1;
    end
end


%% filter the low-D data
for i = 1: numel(dataForJPCA)
    %    sigToFilter = resample(double(dataForJPCA(i).A), binsize_rescaled, 1);
    sigToFilter = dataForJPCA(i).A';
    sigFiltered = bandpassFilter_singleTrial( double(sigToFilter), 14, 9, 100);
    %    dataForJPCA(i).A = resample(sigFiltered', 1, binsize_rescaled);
    dataForJPCA(i).A = sigFiltered';
end


%% perform jPCA
jPCA_params.softenNorm = 5;
jPCA_params.suppressBWrosettes = true;
jPCA_params.suppressHistograms = true;
%%
times = 500:10:730;
jPCA_params.numPCs = 8;
[Projection, Summary] = jPCA(dataForJPCA, times, jPCA_params);

%%
plotParams.planes2plot = [ 1 2 ];
phaseSpace(Projection,Summary, plotParams);

%% plot progression on certain jPCs
splitplot = 0;
clear end
trail_length = 3;
%for nt = 46
for nt = 1: 2: numel( timePoints_jPCA )
    startind = max(1, nt - trail_length );
    close all
    f1 = figure;
    if ~splitplot
        set(gcf,'position',[204    85   924   793]);
    else
        set(gcf, 'position', [4         271        1349         603]);
    end
    for i = 1:numel(Projection)
        if (RT(i) < median(RT) + 0.05 && RT(i) > median(RT) - 0.05) || (DL(i) < 80 || DL(i) > 90)
            %if DL(i) < 80 || DL(i) > 90 || RT(i) > median(RT) + 0.08 || RT(i) < median(RT) -0.08 % this is for near median
            continue;
        end
        
        if RT(i)>=median(RT) + 0.05
            thisColor = [0 0 1];
            plotind = 1;
        elseif RT(i)<=median(RT) - 0.05
            thisColor = [1 0 0];
            plotind = 2;
        end

        if splitplot
            subplot(1,2,plotind);
        end
        
        h = Plot.patchline( Projection(i).projAllTimes( startind : nt, 1), Projection(i).projAllTimes( startind : nt, 2),'edgecolor',thisColor,'linewidth',0.3,'edgealpha',0.5, 'facealpha', ...
                            0.5);
        %set( h, 'edgecolor', linecolor );
        %                set( h, 'edgecolor', linecolor * thisRtRatio );
        %set( h, 'edgecolor', linecolor );
        %set( h, 'facealpha', 0.5, 'edgealpha', 0.5 );
        %set( h, 'linewidth', 0.3);
        hold on;
        h2 = plot(Projection(i).projAllTimes( nt, 1) , Projection(i).projAllTimes( nt, 2), '.');
        set(h2, 'markersize', 5);
        %set(h2, 'color', linecolor );
        %set(h2, 'color', linecolor*thisRtRatio );
        set(h2, 'color', thisColor)
        
            %        set(h2, 'color', 'k' );
        %set( h, 'edgealpha', 0.1 + thisRtRatio/1.2 );
        hold on
    end

    if splitplot
        p=[];
        for nn= 1:2
            p(nn) = subplot(1,2,nn);
        end
    else
        p = gca();
    end

    for np = 1:2
        title(p(np), [int2str(Projection(1).allTimes(nt)) ' ms'])
        xlim(p(np), [-0.12 0.11])
        ylim(p(np), [-0.08 0.2])
    end
    
        %xlim([-0.12 0.11]) % for no filter
        %ylim([-0.08 0.2]) % for no filter
    %xlim([-0.06 0.06])
    %ylim([-0.06 0.06])
    %    title([int2str(Projection(1).allTimes(nt)) ' ms'])
    filename= sprintf( '%s%04g.png', outdir, nt);
    img = getframe(gcf);
    savepng( img.cdata, filename, 4 );
    %h_SFN = gcf();
    %printpdf(h_SFN,int2str(nt) )
    % savefig(h_SFN, int2str(nt));
    pause(0.1);
end

%%

%% plot progression on certain jPCs
clear end
trail_length = 3;
%for nt = 46
for nt = 1: 2: numel( timePoints_jPCA )
    startind = max(1, nt - trail_length );
    close all
    f1 = figure;
    for i = 1:numel(Projection)
        %if RT(i)<=median(RT) + 0.05
        %    continue;
        %end
        
        %        if (RT(i) < median(RT) + 0.05 && RT(i) > median(RT) - 0.05) || (DL(i) < 83 || DL(i) > 85)
        %           continue;
        %        end
        
        h = Plot.patchline( Projection(i).projAllTimes( startind : nt, 1), Projection(i).projAllTimes( startind : nt, 2),'edgecolor',[0.8500 0.5 0.0980],'linewidth',0.3,'edgealpha',0.5, 'facealpha', ...
                            0.5);
        %set( h, 'edgecolor', linecolor );
        %                set( h, 'edgecolor', linecolor * thisRtRatio );
        %set( h, 'edgecolor', linecolor );
        %set( h, 'facealpha', 0.5, 'edgealpha', 0.5 );
        %set( h, 'linewidth', 0.3);
        hold on;
        h2 = plot(Projection(i).projAllTimes( nt, 1) , Projection(i).projAllTimes( nt, 2), '.');
        set(h2, 'markersize', 5);
        %set(h2, 'color', linecolor );
        %set(h2, 'color', linecolor*thisRtRatio );
        %if RT(i)>=median(RT) + 0.05            
        %    set(h2, 'color', 'k' );
        %elseif RT(i)<=median(RT) - 0.05   
        %    set(h2, 'color', 'r' );
        %end
        
        set(h2, 'color', 'k' );
        %set( h, 'edgealpha', 0.1 + thisRtRatio/1.2 );
        hold on
    end
    xlim([-0.12 0.11]) % for no filter
    ylim([-0.08 0.2]) % for no filter
    %xlim([-0.06 0.06])
    %ylim([-0.06 0.06])
    %timeName = (DL(i)+timePoints_jPCA(nt))*10;
    title([int2str(Projection(1).allTimes(nt)) ])
    filename= sprintf( '%s%04g.png', outdir, nt);
    img = getframe(gcf);
    savepng( img.cdata, filename, 4 );
    %h_SFN = gcf();
    %printpdf(h_SFN,int2str(nt) )
    % savefig(h_SFN, int2str(nt));
    pause(0.1);
end


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

numTrials = length(dataForProjection);
allTimes = timePoints_jPCA*binsize_rescaled;
analyzeIndices = ismember(round(allTimes), times);
analyzeMask = repmat(analyzeIndices,numTrials,1);  % used to mask bigA
clear dataForProjection
dataForProjection(totalTrialsToKeep).A = 0;
FRandSpiking = [];
ind = 1;
for nday = 1:numel( alf )
%for nday = [1 3]
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        dataForProjection(ind).A = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        dataForProjection(ind).A = dataForProjection(ind).A';        
        dataForProjection(ind).times = timePoints_jPCA*binsize_rescaled;
        dataForProjection(ind).times = dataForProjection(ind).times';
        ind = ind + 1;
    end
end
%%
clear times
bigA = vertcat(dataForProjection.A);  % append conditions vertically
%bigA = vertcat(dataForJPCA.A);  % append conditions vertically
normFactors = Summary.preprocessing.normFactors;
meanFReachNeuron = Summary.preprocessing.meanFReachNeuron;
bigA = bsxfun(@times, bigA, 1./normFactors);  % normalize
sumA = 0;
for c = 1:numTrials
    sumA = sumA + bsxfun(@times, dataForProjection(c).A, 1./normFactors);  % using the same normalization as above
    %sumA = sumA + bsxfun(@times, dataForJPCA(c).A, 1./normFactors);  % using the same normalization as above
end
meanA = sumA/numTrials;
bigA = bigA-repmat(meanA,numTrials,1);
smallA = bigA(analyzeMask,:);
proj_allData = bsxfun(@minus, smallA, meanFReachNeuron) * Summary.jPCs_highD; % projection of all the data
index1 = 1;
numAnalyzedTimes = size(smallA,1)/numTrials;
for c = 1:numTrials
    index1b = index1 + numAnalyzedTimes -1;  % we will go from index1 to this point
    newData_projection(c).proj = proj_allData(index1:index1b,:);
    newData_projection(c).times = dataForProjection(1).times(analyzeIndices);
    index1 = index1+numAnalyzedTimes;
end
%%
for i = 1:numTrials
    plot(newData_projection(i).proj(:, 5), newData_projection(i).proj(:, 6));
    hold on
end
%%
    


%% project data to jPCs without filtering
clear dataToProject
dataToProject(totalTrialsToKeep).A = 0;
FRandSpiking = [];
allTimes = timePoints_jPCA*binsize_rescaled;
analyzeIndices = ismember(round(allTimes, times);
ind = 1;
for nday = 1:numel( alf )
%for nday = [1 3]
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        bigA = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        dataToProject(ind).A = bigA(analyzeIndices,:)';        
        dataToProject(ind).times = dataForJPCA(ind).times';
        ind = ind + 1;
    end
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

%% center the data
% for LFP
%for i = 1:numel(chopResampled_lfp)
%    meanThisTrial = meanchopResampled_lfp(i).lfps

%% For LFP data
maxLag = 80;
ind = 0;
%crossCorr{ numel( alf ) } = [];
crossCorr = {};
crossCorr_shuffled = {};
for nday = 1:numel( alf )
    trialIds = ( ind + 1 ):( ind + trialsToKeepEachDay( nday ));
    for n = 1 : size(lfp{ nday }(1).lfps, 1)
        for nTrial = 1 : numel( trialIds )
            trialsToChoose = trialIds;
            trialsToChoose(nTrial) = [];
            trialForShuffle = randsample(trialsToChoose, 1);
            for dimId = 1 : 4
                fieldToStore = ['dim' num2str(dimId)];
                sig1 = Projection( trialIds( nTrial ) ).proj( :, dimId)';
                sig2 = chopResampled_lfp( trialIds( nTrial ) ).lfps( n, :);
                sig3 = chopResampled_lfp( trialForShuffle ).lfps( n, :);
                sig1 = sig1-mean(sig1);
                sig2 = sig2-mean(sig2);
                sig3 = sig3-mean(sig3);
                
                crossCorr{ nday }( n ).( fieldToStore )( nTrial, :) = xcorr( sig1, sig2, maxLag/binsize_rescaled);
                crossCorr_shuffled{ nday }( n ).( fieldToStore )( nTrial, :) = xcorr(sig1, sig3, maxLag/binsize_rescaled);
            end
        end
    end
    ind = ind + trialsToKeepEachDay( nday );
end


%% For spiking data
maxLag = 80;
ind = 0;
%crossCorr{ numel( alf ) } = [];
crossCorr = {};
crossCorr_shuffled = {};
for nday = 1:numel( alf )
    trialIds = ( ind + 1 ):( ind + trialsToKeepEachDay( nday ));
    for n = 1 : size(alf{ nday }(1).spikes, 1)
        for nTrial = 1 : numel( trialIds )
            trialsToChoose = trialIds;
            trialsToChoose(nTrial) = [];
            trialForShuffle = randsample(trialsToChoose, 1);
            for dimId = 1 : 4
                fieldToStore = ['dim' num2str(dimId)];
                crossCorr{ nday }( n ).( fieldToStore )( nTrial, :) = xcorr(FRandSpiking( trialIds( nTrial ) ).spiking( n, :), Projection( trialIds( nTrial ) ).proj( :, dimId)', maxLag/binsize_rescaled);
                crossCorr_shuffled{ nday }( n ).( fieldToStore )( nTrial, :) = xcorr( FRandSpiking( trialForShuffle ).spiking( n, :), Projection(trialIds( nTrial ) ).proj( :, dimId)', maxLag/binsize_rescaled);
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
        shadedErrorBar([], crossCorr{ nday }( n ).dim1, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim1, 1)) }, 'lineProps', '-r');
        shadedErrorBar([], crossCorr_shuffled{ nday }( n ).dim1, {@mean, @(x) std(x)./sqrt(size(crossCorr_shuffled{ nday }( n ).dim1, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});        
        set(gca,'XTick',[1 0.5*size(crossCorr{ nday }( n ).dim1, 2) size(crossCorr{ nday }( n ).dim1, 2)]);
        timeStrLeft = ['-', timeLagStr];
        timeStrRight = ['+', timeLagStr];
        set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
        set(gca,'XLim',[0 size(crossCorr{ nday }( n ).dim1, 2)]);
        ylabel('cross-correlation');
        title('Dim 1');
        hold on
        
        sp2 = subplot(2,2,2);
        shadedErrorBar([], crossCorr{ nday }( n ).dim2, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim2, 1)) }, 'lineProps', '-r');
        shadedErrorBar([], crossCorr_shuffled{ nday }( n ).dim2, {@mean, @(x) std(x)./sqrt(size(crossCorr_shuffled{ nday }( n ).dim2, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});        
        set(gca,'XTick',[1 0.5*size(crossCorr{ nday }( n ).dim2, 2) size(crossCorr{ nday }( n ).dim2, 2)]);
        timeStrLeft = ['-', timeLagStr];
        timeStrRight = ['+', timeLagStr];
        set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
        set(gca,'XLim',[0 size(crossCorr{ nday }( n ).dim2, 2)]);
        ylabel('cross-correlation');
        title('Dim 2');
        hold on

        sp3 = subplot(2,2,3)
        shadedErrorBar([], crossCorr{ nday }( n ).dim3, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim3, 1)) }, 'lineProps', '-r');
        shadedErrorBar([], crossCorr_shuffled{ nday }( n ).dim3, {@mean, @(x) std(x)./sqrt(size(crossCorr_shuffled{ nday }( n ).dim3, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});        
        set(gca,'XTick',[1 0.5*size(crossCorr{ nday }( n ).dim3, 2) size(crossCorr{ nday }( n ).dim3, 2)]);
        timeStrLeft = ['-', timeLagStr];
        timeStrRight = ['+', timeLagStr];
        set(gca,'XTickLabels',{timeStrLeft,'0',timeStrRight});
        set(gca,'XLim',[0 size(crossCorr{ nday }( n ).dim3, 2)]);
        ylabel('cross-correlation');
        title('Dim 3');
        hold on

        sp4 = subplot(2,2,4)
        shadedErrorBar([], crossCorr{ nday }( n ).dim4, {@mean, @(x) std(x)./sqrt(size(crossCorr{ nday }( n ).dim4, 1)) }, 'lineProps', '-r');
        shadedErrorBar([], crossCorr_shuffled{ nday }( n ).dim4, {@mean, @(x) std(x)./sqrt(size(crossCorr_shuffled{ nday }( n ).dim4, 1)) }, 'lineProps', {'Color',[0.6, 0.6, 0.6]});
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