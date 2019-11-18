%% make a place to store output videos
outdir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/reRun180614_20190620/lowD_traj/cue123/cond_avg/';
if ~isdir( outdir )
    mkdir( outdir );
end

cd(outdir)

% number of trials for each day
numTrialsTot = cellfun( @numel, alf );


% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial


for nday = 1 : numel( alf )
    isCorrectArray{ nday } = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nday }.arrayShapesCorrect);
    %isCueLoc1{ nday } = UEs{ nday }.cueLoc == 1;
    trialsToKeep{ nday } = isCorrectArray{ nday }; % & isCueLoc1{ nday };
    cueLocs{ nday } = unique(UEs{ nday }.cueLoc);

    for nc = 1 : numel( cueLocs{ nday } )
        trialsByCueLoc{ nday }{nc} = find( trialsToKeep{ nday } & (UEs{ nday }.cueLoc==cueLocs{ nday }(nc)));
        rtsByCueLoc{ nday }{nc} = UEs{ nday }.rt( trialsByCueLoc{ nday }{nc} );    
    end
end

%

% concatenate all the factors

% % this window was used for finding oscillations during arrayDelay in
% %         factors 6 7 8 for session 6
%window = round( [-300  00] / binsize_rescaled );

window = round([0 400]/binsize_rescaled);

%window = round( [-1200 : 500] / binsize );

%whichfieldDimred = 'arrayDim';
%whichfieldDimred = 'arrayOnset';
whichfieldDimred = 'cueOnset';

%whichfieldPlot = 'arrayDim';
%newWindow = round( [-900  650] / binsize_rescaled );


%whichfieldPlot = 'arrayOnset';
%newWindow = round( [0  1200] / binsize_rescaled );


whichfieldPlot = 'cueOnset';
newWindow = round( [0  900] / binsize_rescaled );

%% condition averaging and then PCA
timePoints = window(1):window(2);
numBins = numel( timePoints );
numFactors = size( alf{ nday }(1).rates, 1);
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );
allFactors_1 = zeros( numFactors, numBins);
allFactors_3 = zeros( numFactors, numBins);

%
c1_count = 0;
c3_count = 0;
for nday = 1 : numel( alf)
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        if UEs{nday}.cueLoc( ntr ) == 1
            allFactors_1 = allFactors_1 + alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
            c1_count = c1_count + 1;
        else
            allFactors_3 = allFactors_3 + alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
            c3_count = c3_count + 1;
        end
    end
end

allFactors_1 = allFactors_1/c1_count;
allFactors_3 = allFactors_3/c3_count;
allFactors = [allFactors_1, allFactors_3];

%% do PCA on single-trials

timePoints = window(1):window(2);
numBins = numel( timePoints );
numFactors = size( alf{ nday }(1).rates, 1);
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );
allFactors = zeros( numFactors, numBins * totalTrialsToKeep );

%
ind = 1;
for nday = 1 : numel( alf)
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        allFactors( :, (0:numBins-1) + ind ) = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        ind = ind + numBins;
    end
end

%% do pca
meanFactors = mean( allFactors' );

[pca_proj_mat, pc_data] = pca( allFactors', 'NumComponents', 30);



%% %  % plot
cmap = lines();
% which pcs to plot?
p2p = [1 2 3];
%paxis=[-2.1342   -0.6698   -2.5038    1.0200   -1.3868    1.3984];
%paxis = [-1.7655   -0.7427   -1.8530    0.4962   -0.5494    1.3074];

%paxis([1 3 5] ) = paxis([1 3 5])*2;
%paxis = [-1.4930   -0.6525   -1.7117    1.1907   -1.0502    2.2128];

%p2p = [4 5 6];
%paxis = [-1.2255    0.8949   -1.0455    0.9195   -1.1087    1.3748];



%p2p = [6 7 8];
p2p = [1 2 3];
pmin1 = min(pc_data(:, p2p(1)));
pmax1 = max(pc_data(:, p2p(1)));
pmin2 = min(pc_data(:, p2p(2)));
pmax2 = max(pc_data(:, p2p(2)));
pmin3 = min(pc_data(:, p2p(3)));
pmax3 = max(pc_data(:, p2p(3)));
paxis = [pmin1 pmax1 pmin2 pmax2 pmin3 pmax3];

% p2p = [10 11 12 ];
%p2p = [13 14 15 ];


window = newWindow(1):newWindow(2);
numBins = numel( window );

plotFactors_1 = zeros( numFactors, numBins);
plotFactors_3 = zeros( numFactors, numBins);

%
c1_count = 0;
c3_count = 0;
for nday = 1 : numel( alf)
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        if UEs{nday}.cueLoc( ntr ) == 1
            plotFactors_1 = plotFactors_1 + alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldPlot ) + window );
            c1_count = c1_count + 1;
        else
            plotFactors_3 = plotFactors_3 + alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldPlot ) + window );
            c3_count = c3_count + 1;
        end
    end
end

plotFactors_1 = plotFactors_1/c1_count;
plotFactors_3 = plotFactors_3/c3_count;


rtRatios = [];

splitplot = 0;
%nplots=2; % cp
clf;
if ~splitplot
    set(gcf,'position',[204    85   924   793]);
else
    set(gcf, 'position', [4         271        1349         603]);
end


% trail length in bins
%trail_length = 20;% cp
trail_length = 20;

%%
%
for nt = 1 : 2 :numel( window )
    %for nt = 26
    %for nt = 1:5:80
    clf;
    for icond = 1:2
        if icond ==1
            frep = plotFactors_1;
            cueInd = 1;
        else
            frep = plotFactors_3;
            cueInd = 3;
        end
        
        if splitplot
            if icond == 1
                plotind = 1;
            else
                plotind = 2;
            end
        end
        
        % mean center
        frep = frep - repmat( meanFactors(:), 1, numBins );

        % project this data
        dim_reduced_data = pca_proj_mat' * frep;

        linecolor = cmap( cueInd, : );

        plot_mean = 0;

        if ~plot_mean
            % plot the trace

            startind = max(1, nt - trail_length );

            if splitplot
                subplot(1,nplots,plotind);
            end
            h = Plot.patchline( dim_reduced_data( p2p(1), startind : nt), ...
                                dim_reduced_data( p2p(2), startind : nt), ...
                                dim_reduced_data( p2p(3), startind : nt) );
            %set( h, 'edgecolor', linecolor );
            %                set( h, 'edgecolor', linecolor * thisRtRatio );
            set( h, 'edgecolor', linecolor );
            set( h, 'facealpha', 0.5, 'edgealpha', 0.5 );
            hold on;


            h2 = plot3( dim_reduced_data( p2p(1), nt ), ...
                        dim_reduced_data( p2p(2), nt ), ...
                        dim_reduced_data( p2p(3), nt ), '.');
            set(h2, 'markersize', 10);
            %set(h2, 'color', linecolor );
            %set(h2, 'color', linecolor*thisRtRatio );
            set(h2, 'color', linecolor );
            %set( h, 'edgealpha', 0.1 + thisRtRatio/1.2 );
        else
            % plot the mean
            drdmean = mean( dim_reduced_data' );
            h = plot3( drdmean( p2p(1) ), ...
                       drdmean( p2p(2) ), ...
                       drdmean( p2p(3) ), 'o' );
            set( h, 'markerfacecolor', linecolor, 'markeredgecolor', linecolor );
            hold on;
        end
    end
    
    if splitplot
        p=[];
        for nn= 1:nplots
            p(nn) = subplot(1,nplots,nn);
        end
    else
        p = gca();
    end

    clear l l2
    for np = 1:numel(p)
        title(p(np), window( nt ) * binsize_rescaled );
        xlabel(p(np), p2p(1) );
        ylabel(p(np), p2p(2) );
        zlabel(p(np), p2p(3) );
        % this works for dim [6 7 8];
        %set(p(np), 'view', [-23.6000    7.6000]);
        
        %set(p(np), 'view', [ -198.0000   -2.8000]);
        %set(p(np), 'view', [-30.4000 7.6000]);
        %set(p(np), 'view', [-6.8000 74.8000]);
        %set(p(np), 'view', [-6.8000 90]);
        set(p(np), 'view', [-39.2000   19.6000]); % many videos made using this view
                                                  %axis(p(np), [ -1.36    0.4   -0.79    0.63   -1.3615    0.9]);
                                                  %axis(p(np), [-1.5    0.4   -1.10    1.7203   -1.9    0.6567]);% this works well for cueOnset
                                                  %axis(p(np), [-2.5    1   -1.5   2   -2    2]); % SFN cueOnset
                                                  %axis(p(np), [ -2    2   -2   2   -2    2]); % SFN arrayOnset
        axis(p(np), [-5    5   -5   5   -5    5]); % SFN cueOnset
        
        %if ~isempty( paxis )
        %    axis(p(np), paxis );
        %end
        %else
        %    axis('tight');
        %end
        set(p(np),'LooseInset',get(p(np),'TightInset'));
        set(p(np), 'Position', get(p(np), 'OuterPosition') - ...
                   get(p(np), 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        tmp = axis(p(np));
        l(np, :) = tmp;
    end

    if splitplot
        l2([1 3 5]) = min( l(:, [1 3 5] ) );
        l2([2 4 6]) = max( l(:, [2 4 6] ) );
        axis(p, l2 );
    end
    filename= sprintf( '%s%04g.png', outdir, nt);
    img = getframe(gcf);
    savepng( img.cdata, filename, 4 );
    %h_SFN = gcf();
    %printpdf(h_SFN,int2str(nt) )
    % savefig(h_SFN, int2str(nt));
    pause(0.1);
end     