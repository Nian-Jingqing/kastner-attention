function [ dataBinned ] = rebinChannel( data, origSR, newSR )
% warn about upsampling
    flipFlag = false;
    if origSR < newSR
        disp( 'WARNING: Upsampling data ...' )
    end
    % put data in row matrix
    if size( data, 2 ) < size( data, 1 )
        data = data';
        flipFlag = true;
    end

    % find bin size in ms
    binSize = 1000 / newSR;
    % trim any data that won't fit in an even number of bins
    dataDim = size( data, 1 );
    % find dataLength in ms
    dataLength = ( size( data, 2 ) / origSR ) * 1000;
    numBins = floor( dataLength / binSize );
    pointsToKeep = numBins * binSize;
    data = data( :, 1:pointsToKeep );

    
    % reshape, sum and squeeze
    d = reshape( full( data ), dataDim, binSize, numBins );
    dataBinned = squeeze( sum( d, 2 ) );
    
    % if the original data is 1-dimensional, there are some problems with squeeze
    %  correction for that is below
    if dataDim == 1
        dataBinned = dataBinned(:)';
    end

    if flipFlag
        dataBinned = dataBinned';
    end
    
end
