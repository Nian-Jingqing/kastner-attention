%% set path
outdir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/reRun180614_20190620/jPCA/projected/tmp/';
if ~isdir( outdir )
    mkdir( outdir );
end
cd(outdir)
%% pre-processing
% number of trials for each day
numTrialsTot = cellfun( @numel, alf );

%% selecting certain condition for analysis
for nday = 1 : numel( alf )
    isCorrectArray{ nday } = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nday }.arrayShapesCorrect);
    isCueLoc3{ nday } = UEs{ nday }.cueLoc == 3;
    trialsToKeep{ nday } = isCorrectArray{ nday } & isCueLoc3{ nday };
    cueLocs{ nday } = unique(UEs{ nday }.cueLoc);
    for nc = 1 : numel( cueLocs{ nday } )
        trialsByCueLoc{ nday }{nc} = find( trialsToKeep{ nday } & (UEs{ nday }.cueLoc==cueLocs{ nday }(nc)));
        rtsByCueLoc{ nday }{nc} = UEs{ nday }.rt( trialsByCueLoc{ nday }{nc} );    
    end
end

%% applying highpass filter
alf = ala;

for nday = 1 : numel( alf)
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        alf{ nday }( ntr ).rates(:, alf{ nday }( ntr ).cueOnset : alf{ nday }( ntr ).arrayDim+280/8) = highpassFilter_singleTrial( double(alf{ nday }( ntr ).rates(:, alf{ nday }( ntr ).cueOnset : alf{ nday }( ntr ).arrayDim+280/8)), 3, 125);
        %tmp  = highpassFilter_singleTrial( double(alf{ nday }( ntr ).rates(:, alf{ nday }( ntr ).arrayOnset : alf{ nday }( ntr ).arrayDim+280/8)), 3, 125);
        %alf{ nday }( ntr ).rates(:, alf{ nday }( ntr ).arrayOnset : alf{ nday }( ntr ).arrayDim+280/8) = normalize(tmp, 'noCenter');
    end
end


%% selecting certain time window for analysis
%window = round( [-300  0] / binsize_rescaled );
%whichfieldDimred = 'arrayDim';

% selecting the window for doing jPCA
whichfieldJPCA = 'arrayDim';
window_jPCA = round( [-600  200] / binsize_rescaled );

% selecting the window for subsequent plotting
whichfieldPlot = 'arrayOnset';
newWindow = round( [-200  1200] / binsize_rescaled );

%% gather relevant data information
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );
for nday = 1 : numel( alf )
    nNeurons(nday) = size(alf{nday}(1).FR, 1);
end
numFactors = size( alf{ 1 }(1).rates, 1);
timePoints_jPCA = window_jPCA(1):window_jPCA(2);

%% perform jPCA
clear dataForJPCA
dataForJPCA(totalTrialsToKeep).A = 0;
% FRandSpiking = [];
ind = 1;
RT = zeros(1, totalTrialsToKeep);
DL = zeros(1, totalTrialsToKeep);
for nday = 1:numel( alf )
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        dataForJPCA(ind).A = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldJPCA ) + timePoints_jPCA );
        dataForJPCA(ind).A = dataForJPCA(ind).A';        
        dataForJPCA(ind).times = timePoints_jPCA*binsize_rescaled;
        dataForJPCA(ind).times = dataForJPCA(ind).times';
        RT(ind) = alf{ nday }( ntr ).rt;
        DL(ind) = alf{ nday }( ntr ).arrayDim - alf{ nday }( ntr ).arrayOnset;
        ind = ind + 1;
    end
end


%% filter the low-D data
for i = 1: numel(dataForJPCA)
    sigToFilter = dataForJPCA(i).A';
    sigFiltered = bandpassFilter_singleTrial( double(sigToFilter), 8, 3, 125);
    dataForJPCA(i).A = sigFiltered';
end

%% perform jPCA
jPCA_params.softenNorm = 5;
jPCA_params.suppressBWrosettes = true;
jPCA_params.suppressHistograms = true;
%%
times = -296:8:-32;
jPCA_params.numPCs = 8;
[Projection, Summary] = jPCA(dataForJPCA, times, jPCA_params);

%%
plotParams.planes2plot = [ 1 2 ];
phaseSpace(Projection,Summary, plotParams);

%% plot progression on certain jPCs
cmap = lines();
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
        %if (RT(i) < median(RT) + 0.05 && RT(i) > median(RT) - 0.05) || (DL(i) < 80 || DL(i) > 90)
            %if DL(i) < 80 || DL(i) > 90 || RT(i) > median(RT) + 0.08 || RT(i) < median(RT) -0.08 % this is for near median
            %   continue;
            %end
        
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

        linecolor = cmap(2,:);
        h = Plot.patchline( Projection(i).projAllTimes( startind : nt, 1), Projection(i).projAllTimes( startind : nt, 2),'edgecolor',linecolor,'linewidth',0.3,'edgealpha',0.5, 'facealpha', ...
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
        set(h2, 'color', 'k')
        
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

    for np = 1
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