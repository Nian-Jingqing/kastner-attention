function [ ] = plotNeural( run, par, plotTime, dataSize, trialNum, writeFlag, outputDir, filePrefix, trialId )
    pca = run.Data.PCA( :, plotTime );
    spikes = run.Data.Spikes( :, plotTime );
    smoothedSpikes = run.Data.SmoothedSpikes( :, plotTime );
    factors = run.Data.Factors( :, plotTime );
    nFactors = size( factors, 1 );

    rates = run.Data.Rates( :, plotTime );
    nRates = size( rates, 1 );
    sizeRatio = nFactors / nRates;
    sizeRatio = 0.4;
    vertSpace = 0.0;
    figure()

    % plot binned spikes
    subaxis( 5, 1, 1, 'SpacingVert', vertSpace)
    hSpike = imagesc( spikes );
    box('off')
    %ylabel( 'Bin. Spikes')
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
    spikePlotPos = hSpike.Parent.Position;
    
    % plot smoothed spikes
    subaxis( 5, 1, 2, 'SpacingVert', vertSpace )
    hSmooth = imagesc( smoothedSpikes );
    %ylabel( 'Smooth. Spikes')
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
    smoothPlotPos = hSmooth.Parent.Position;
    
    subaxis( 5, 1, 3, 'SpacingVert', vertSpace )
    hRate = imagesc( rates );
    %ylabel( 'Rates')
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
    ratePlotPos = hRate.Parent.Position;
    if writeFlag
        filename = strcat( filePrefix, '-trial-', num2str( trialId ), '-', 'neural' );
        filepath = fullfile( outputDir, filename );
        print( filepath,'-dpdf','-r300', '-bestfit' )
        close all
    end
    %{
    subaxis( 5, 1, 4, 'SpacingVert', vertSpace )
    hPCA = imagesc( pca );
    %ylabel( 'PCA')
    currentPos = Cherian.utils.changeSubaxisPlotSize( hPCA, ratePlotPos, vertSpace, sizeRatio );
    hPCA.Parent.Position = currentPos;
    pcaPlotPos = hPCA.Parent.Position;
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )

    subaxis( 5, 1, 5, 'SpacingVert', vertSpace )
    hFac = imagesc( factors );
    %ylabel( 'Factors');
    currentPos = Cherian.utils.changeSubaxisPlotSize( hFac, pcaPlotPos, vertSpace, sizeRatio );
    hFac.Parent.Position = currentPos;
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
    %}
    

    %{
    nTicks = [ 0:100:round( dataSize, -2 ) ];
    nTicks = nTicks + 1;
    set( gca, 'XTick', nTicks );
    nTicks = nTicks - 1;
    for j = 1:numel(nTicks)
        XTickLabels{ j } = num2str( ( nTicks( j ) * par.spikeBinMs ) );
    end
    set( gca, 'XTicklabel', XTickLabels )
    %}
end
