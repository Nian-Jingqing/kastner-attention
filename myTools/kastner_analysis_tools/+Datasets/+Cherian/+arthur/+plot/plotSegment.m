function [ ] = plotSegment( data, sampleRate, segSizeMs, chanIdx, segIdx )
    % get number of sample in segSize
    segSize = segSizeMs / ( 1000 / sampleRate );
    % extract a channel to get some info to calculate seg size
    channel  = data( chanIdx, : );
    % determine length of channels
    channelLength = size( channel, 2 );
    % calculate seg size
    nSegs = round( channelLength / segSize );
    % extract segment
    segment = channel( ( segIdx-1 )*segSize+1 : segIdx*segSize );

    % initialize plot
    Plot.blankFigure();
    opacity = 0.7;
    x = 1 : size( segment, 2 );
    y = segment;
    smoothColor = [ 83 229 146 ]./ 255; % wintergreen
    h = Plot.patchline( x, y, ...
                        'LineWidth', 2, ...
                        'edgeAlpha', opacity, ...
                        'edgeColor', smoothColor, ...
                        'AlignVertexCenters', 'off' );
    figSizeX = 10;
    figSizeY = 5;
    figure_prop_name = { 'PaperPositionMode','units','Position' };
    figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };
    set( gcf, figure_prop_name, figure_prop_val )
end
