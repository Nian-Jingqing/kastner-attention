function [ ] = plotCharacteristics( data, segLength, yLabel, plotChar, outputDir, filePrefix, writeFlag )
    figure();
    cmap = hsv;
    nMuscles = size( data, 1 );
    % find number of segments
    nSegs = 1:size( data( 1, : ), 2 );
    opacity = 0.7;
    for iMuscle = 1 : nMuscles
        % segment characteristic for each muscle
        iMuscleChar = data( iMuscle , : );
        if iMuscle == 5
            keyboard
        end
        
        clrind = floor( (iMuscle - 1) / nMuscles * size(cmap, 1) ) + 1;
        h = Plot.patchline( nSegs, iMuscleChar, ...
                            'LineWidth', 2, ...
                            'edgeAlpha', opacity, ...
                            'AlignVertexCenters', 'off' );
        legend( 'cP', 'lD', 'aD', 'mD', 'mB','lB','lT','bR','eR','fR' )
        hold on;
        set( h, 'edgecolor', cmap(clrind, :) );
        xlabel( 'Seg Number' )
        ylabel( yLabel )
        %h = plot( nSegs, iMuscleChar, 'o', 'MarkerSize', 5 );
        %set( h, 'MarkerEdgeColor', cmap(clrind, :) );
        %set( h, 'MarkerFaceColor', cmap(clrind, :) );
        figSizeX = 10;
        figSizeY = 5;
        figure_prop_name = { 'PaperPositionMode','units','Position' };
        figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };
        plotTitle = ['Segment Length ', plotChar, ':', num2str( segLength ),'s' ];
        outTitle = [ plotChar, '-segSize-', num2str( segLength ), 's' ];
        title( plotTitle )
        set( gcf, figure_prop_name, figure_prop_val )
        orient( gcf, 'landscape' )
    end
    hold off
    if writeFlag
        outTitle = lower( outTitle );
        fileSpec = strcat( '-', outTitle );
        filename = strcat( filePrefix, fileSpec );
        filepath = fullfile( outputDir, filename );
        print( filepath,'-dpdf','-r300','-bestfit' )
        print( filepath,'-dpng','-r300' )
        close all
    end

end
