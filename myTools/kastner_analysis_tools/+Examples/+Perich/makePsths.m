%%
ddir = '/snel/share/share/data/Miller_Perich/';
matname = 'Chewie_20161011_CO_FF_BL_001_stripped.mat';
%matname = 'Chewie_20150319_CO_CS_BL_001_stripped.mat';

fname = fullfile( ddir, matname );

%%
[C startInds stopInds trialstruct ] = Datasets.Perich.loadAndProcess( fname );

%%
% smooth the spiketrains with a gaussian kernel
C.smoothField( 'spikes', 'y_smoothed', 100 );

%%
clear r
% turn into a trialized (R) struct
r = R.Rstruct( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

%% trim to successful trials only
outcomes = [r.r.result];
r.r = r.r( outcomes == 'R' );



%% get some individual trials
alignedTarg = r.getAligned('conditionID', {'y_smoothed', 'spikes', 'kin'}, [ 50 500 ], 'targetOnset' );
alignedGo = r.getAligned('conditionID', {'y_smoothed', 'spikes', 'kin'}, [ 200 800 ], 'goCue' );


%% get a colormap
cmap = hsv;
numConds = numel( alignedGo );


%% get some trial-averaged data
targPrePost = [ 50 700 ];
goPrePost = [ 200 800 ];
avgTarg = r.alignAndAverage('conditionID', {'y_smoothed', 'spikes', 'kin'}, targPrePost, 'targetOnset' );
avgGo = r.alignAndAverage('conditionID', {'y_smoothed', 'spikes', 'kin'}, goPrePost, 'goCue' );



%% plot PSTHs for a few neurons
dfield = 'y_smoothed';

t1 = 1 + (-targPrePost(1) : targPrePost(2)-1);
offset = 100;
t2 = 1 + (-goPrePost(1) : goPrePost(2)-1) + max(t1) + ( offset + goPrePost(1) );

cmap = hsv;
numConds = numel( avgGo );
opacity = 0.5;

for icell = 1:100
    Plot.blankFigure( icell );
    for icond = 1:numel( avgTarg )
        clrind = floor( (icond - 1) / numConds * size(cmap, 1) ) + 1;
        h = Plot.patchline( t1, avgTarg( icond ).(dfield)( icell, : ), ...
                            'edgealpha', opacity );
        set( h, 'edgecolor', cmap(clrind, :), 'linewidth', 2 );
        hold on;
        h = Plot.patchline( t2, avgGo( icond ).(dfield)( icell, : ), ...
                            'edgealpha', opacity );
        set( h, 'edgecolor', cmap(clrind, :), 'linewidth', 2  );
    end
    axis('tight');
    keyboard
end


%% plot individual trial kinematics
Plot.blankFigure(101); clf;
opacity = 0.5;
for icond = 1:numel( alignedGo )
    X = alignedGo( icond ).kin;
    clrind = floor( (icond - 1) / numConds * size(cmap, 1) ) + 1;
    for nt = 1:size( X, 1)
        x = squeeze( X( nt, 1, : ) );
        y = squeeze( X( nt, 2, : ) );
        h = Plot.patchline( x,  y, 'edgealpha', opacity );
        hold on;
        set( h, 'edgecolor', cmap(clrind, :) );
    end % nt
    keyboard
end % nc




%% plot a few velocity-based trajectories
Plot.blankFigure(200);
for nc = 1:numel( avgGo )
    pos(nc).X = cumsum( avgGo(nc).kin' );
    plot( pos(nc).X(:, 3), ...
          pos(nc).X(:, 4) );
    hold on;
end

%% plot a few position-based trajectories
Plot.blankFigure(201);
for nc = 1:numel( avgGo )
    pos(nc).X = avgGo(nc).kin;
    plot( pos(nc).X(1, :), ...
          pos(nc).X(2, :) );
    hold on;
end
