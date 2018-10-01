use_rates = false;


% % olapChopped.r.r contains receptive field locations

% % assembled_lfads.r has all the LFADS results
if use_rates
    r_lfads2 = olapChopped.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms, 'rates');
else
    r_lfads2 = olapChopped.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms, 'factors');
end

assembled_lfads2 = R.Rstruct(r_lfads2);
alf = assembled_lfads2.r;

binsize = par.spikeBinMs;

% 'dimStart' is only defined for hold trials
% make a version that's defined for all trials
dimStartAll = nan( size( cueStart ) );
dimStartAll( UE.isHoldTrial ) = dimStart;

rfactor = 2;
binsize = binsize * rfactor;

for ntr = 1:numel( alf )
    if use_rates
        alf( ntr ).rates = log( alf( ntr ).rates );
    end
    rdim = size( alf(ntr).rates );
    tmp = alf(ntr).rates(:, 1:floor( rdim( 2 ) / rfactor ) * rfactor );
    tmp = reshape( tmp, rdim(1), rfactor, floor( rdim( 2 ) / rfactor ) );
    tmp = squeeze( mean( tmp, 2 ) );
    alf(ntr).rates = tmp;
    
    %alf( ntr ).arrayOnsetMS = UE2.arrayOnset( ntr );
    %alf( ntr ).arrayOnset = ceil(UE2.arrayOnset( ntr ) / binsize);
    alf( ntr ).cueOnset = round( cueStart( ntr ) / binsize);
    alf( ntr ).arrayOnset = round( arrayStart( ntr ) / binsize);
    alf( ntr ).arrayDim = round( dimStartAll( ntr ) / binsize);
end

%%

%whichfieldDimred = 'arrayOnset';
%whichfieldPlot = 'arrayOnset';

whichfieldDimred = 'arrayDim';
whichfieldPlot = 'arrayDim';

% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial

isCorrectArray = arrayfun(@(x) strcmp(x, 'HRHR'), UE2.arrayShapesCorrect);
trialsToKeep = isCorrectArray;% & UE2.isHoldTrial;

cueLocs = unique(UE2.cueLoc);

for nc = 1:numel(cueLocs)
    trialsByCueLoc{nc} = find( trialsToKeep & (UE2.cueLoc==cueLocs(nc)));
    rtsByCueLoc{nc} = UE2.rt( trialsByCueLoc{nc} );    
end

% concatenate all the factors

% % this window was used for finding oscillations during arrayDelay in factors 6 7 8
window = round( [-300 : 00] / binsize );

%window = round( [-300 : 400] / binsize );


numBins = numel( window );
numFactors = size( alf(1).rates, 1);
allFactors = zeros( numFactors, numBins * sum( trialsToKeep ) );

ind = 1;
trialsToKeepInds = find( trialsToKeep );
for itr = 1:numel( trialsToKeepInds )
    ntr = trialsToKeepInds( itr );
    allFactors( :, (0:numBins-1) + ind ) = alf( ntr ).rates( :, alf( ntr ).( whichfieldDimred ) + window );
    ind = ind + numBins;
end


% do pca
meanFactors = mean( allFactors' );

[pca_proj_mat, pc_data] = pca( allFactors', 'NumComponents', 15);



%% %  % plot
cmap = lines();
% which pcs to plot?
%p2p = [1 2 3];
%p2p = [4 5 6];
%p2p = [7 8 9];
p2p = [6 7 8];
% p2p = [10 11 12 ];
%p2p = [13 14 15 ];

newWindow = round( [-400  100] / binsize );

window = newWindow(1):newWindow(2);
numBins = numel( window );

arrayOnsetBin = find( window == min( [ window(end) 0 ] ) );

rtRatios = [];

outdir = '/snel/home/cpandar/tmp/kastnervid/';
if ~isdir( outdir )
    mkdir( outdir );
end

splitplot = 0;
clf;
if ~splitplot
    set(gcf,'position',[204    85   924   793]);
else
    set(gcf, 'position', [4         271        1349         603]);
end

for nt = numel(window) %1 : 10 : numel( window )

    clf;
    
    for itr = 1:numel( trialsToKeepInds )
        ntr = trialsToKeepInds( itr );
        cueInd = find( cueLocs == UE2.cueLoc( ntr ) );
        allRtsThisLoc = rtsByCueLoc{ cueInd };
        thisTrialRt = UE2.rt( ntr );

        %if cueInd == 2
        %    continue;
        %end
        delayTime = dimStartAll - arrayStart;
        if delayTime( ntr ) <= 800
            continue
        end

        sortedRts = sort( allRtsThisLoc, 'ascend' );
        thisRtInd = max( find( sortedRts == thisTrialRt ) );
        thisRtRatio = thisRtInd / numel( sortedRts );
        rtRatios(itr) = thisRtRatio;

        % cut out medium RT trials
        if thisRtRatio > 0.25 & thisRtRatio < 0.75
            continue
        end
        %if thisRtRatio > 0.4
        %    continue
        %end

        if splitplot
            if thisRtRatio <= 0.4
                plotind = 1;
            else
                plotind = 2;
            end
        end
        
        % extract factors for this window
        frep = alf( ntr ).rates( :, alf( ntr ).( whichfieldPlot ) + window );
        % mean center
        frep = frep - repmat( meanFactors(:), 1, numBins );

        % project this data
        dim_reduced_data = pca_proj_mat' * frep;

        linecolor = cmap( cueInd, : );

        plot_mean = 0;

        if ~plot_mean
            % plot the trace

            startind = max(1, nt-600 );

            if splitplot
                subplot(1,2,plotind);
            end
            h = Plot.patchline( dim_reduced_data( p2p(1), startind : nt), ...
                                dim_reduced_data( p2p(2), startind : nt), ...
                                dim_reduced_data( p2p(3), startind : nt) );
            %set( h, 'edgecolor', linecolor );
            set( h, 'edgecolor', linecolor * thisRtRatio );
            set( h, 'facealpha', 0.5, 'edgealpha', 0.5 );
            hold on;


            h2 = plot3( dim_reduced_data( p2p(1), nt ), ...
                        dim_reduced_data( p2p(2), nt ), ...
                        dim_reduced_data( p2p(3), nt ), '.');
            set(h2, 'markersize', 20);
            %set(h2, 'color', linecolor );
            set(h2, 'color', linecolor*thisRtRatio );
            
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
        p(1) = subplot(1,2,1);
        p(2) = subplot(1,2,2);
    else
        p = gca();
    end

    clear l l2
    for np = 1:numel(p)
        title(p(np), window( nt ) * binsize );
        xlabel(p(np), p2p(1) );
        ylabel(p(np), p2p(2) );
        zlabel(p(np), p2p(3) );
        set(p(np), 'view', [-39.2000   19.6000]);
        % this works for dim [6 7 8];
        %axis(p(np), [ -1.31    0.4   -0.79    0.5689   -1.3615    0.55]);
        axis('tight');
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
    %savepng( img.cdata, filename, 4 );
    pause(0.1);
end