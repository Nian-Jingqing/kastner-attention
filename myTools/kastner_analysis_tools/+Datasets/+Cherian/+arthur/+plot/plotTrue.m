function [] = plotTrue( run, plotTime, type, targPosition, dataSize, step )
    addpath( '/snel/home/lwimala/bin/analysis_tools/' )
    
    [ x, y ] = Cherian.utils.determineTrueType( run, type, plotTime );
    axes('Units', 'normalized', 'Position', [0 0 1 1])
    if strcmp( type, 'pos' )
    elseif strcmp( type, 'vel' ) == false
        set( gca, 'Xtick', [ ] )
        ylim( [ -20, round( max( y ) ) + 50 ] )
        xlim( [ 0, dataSize ] )
    else
        set( gca, 'Xtick', [ ] )
        xlim( [ 0, dataSize ] )
    end
    
    if ~strcmp( type, 'pos' )
        Plot.patchline( x, y, ...
                        'LineWidth', 4, ...
                        'edgeColor', [ 0 0 0 ], ...
                        'faceAlpha', 1.0, ...
                        'MarkerSize', 8, ...
                        'MarkerFaceColor', 'k', ...
                        'MarkerEdgeColor', 'k', ...
                        'AlignVertexCenters', 'off' )
    else
        x = x( 1, : );
        y = y( 1, : );
        u = x( 2, : );
        v = y( 2, : );

        ylim( [ -15 15 ] )
        xlim( [ -15 15 ] )
        
        tStartColor = [ 28 158 18 ]; % green
        targColor1 = [ 228 191 46 ]; % yellow-gold
        targColor2 = [ 225 105 20 ]; % orange
        targColor3 = [ 171 11 11 ];  % deep red
        colors = { targColor1, targColor2, targColor3 };
        
        Cherian.plot.plotTarget( x( 1 ), y( 1 ), [ 28 158 18 ] );
        
        for i = 1:numel( targPosition )
            targX = targPosition{ i }( 1 );
            targY = targPosition{ i }( 2 );
            Cherian.plot.plotTarget( targX, targY, colors{ i } );
        end
        
        hold on
        quiver( x( 1:step:end ), y( 1:step:end ), u( 1:step:end ), v( 1:step:end ), '-k', ...
                'LineWidth', 4, ...
                'ShowArrowHead', 'on', ...
                'MaxHeadSize', 0.3, ...
                'AutoScaleFactor', 0.5 )
        Plot.patchline( x, y, ...
                        'LineWidth', 6, ...
                        'edgeColor', [ 0 0 0 ], ...
                        'edgeAlpha', 0.45, ...
                        'faceAlpha', 0.0, ...
                        'MarkerSize', 8, ...
                        'MarkerFaceColor', 'k', ...
                        'MarkerEdgeColor', 'k', ...
                        'AlignVertexCenters', 'off' )
    end
    

    set( gca, 'visible', 'off' )    
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
end
