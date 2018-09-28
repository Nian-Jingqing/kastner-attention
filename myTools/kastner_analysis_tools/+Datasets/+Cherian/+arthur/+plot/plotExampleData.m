function [] = plotExampleData( Rsucc, run, par, trialNum)
% add +Plot class from analysis_tools
    addpath( '/snel/home/lwimala/bin/analysis_tools/' )
    
    % extract trial start index
    tSidx = Rsucc( trialNum ).tSidx;
    % extract trial end index
    tEidx = Rsucc( trialNum ).tEidx;

    % extract original data sample rate
    origSR = run.Data.orig_sR;
    % extract new sample rate
    newSR = run.Data.new_sR;

    % covert sample index to new bin width
    tSidx = Cherian.utils.convertSampleIndex( tSidx, origSR, newSR );
    tEidx = Cherian.utils.convertSampleIndex( tEidx, origSR, newSR );
    
    % determine plot time
    plotTime = tSidx:tEidx;
    dataSize = numel( plotTime );

    targPosition = Cherian.utils.findTargetPosition( Rsucc, trialNum );

    step = 20;
    close all;
    %% neural

    % extract binned spikes
    binnedSpikes = run.Data.Spikes( :, plotTime );
    % plot binned spikes
    Plot.blankFigure();
    imagesc( binnedSpikes )
    figSizeX = 10;
    figSizeY = 3;
    figure_prop_name = { 'PaperPositionMode','units','Position' };
    figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };
    set( gcf, figure_prop_name, figure_prop_val )

    %% position

    % plot position
    Plot.blankFigure();
    figSizeX = 10;
    figSizeY = 10;
    figure_prop_name = { 'PaperPositionMode','units','Position' };
    figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };
    set( gcf, figure_prop_name, figure_prop_val )
    
    Cherian.plot.plotTrue( run, plotTime, 'pos', targPosition, dataSize, step );

    %% velocity

    Plot.blankFigure();
    figSizeX = 10;
    figSizeY = 3;
    figure_prop_name = { 'PaperPositionMode','units','Position' };
    figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };
    set( gcf, figure_prop_name, figure_prop_val )

    Cherian.plot.plotVelocity( run, plotTime )

    % emg

    Plot.blankFigure();
    figSizeX = 10;
    figSizeY = 3;
    figure_prop_name = { 'PaperPositionMode','units','Position' };
    figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };
    set( gcf, figure_prop_name, figure_prop_val )

    Cherian.plot.plotEMG( run, plotTime )
    
end
