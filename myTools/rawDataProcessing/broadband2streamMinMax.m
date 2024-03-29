function spikeband = broadband2streamMinMax(inputfile, outputfile, cutOffFreqs,...
                                            options)
    % BROADBAND2STREAMMINMAX
    % 
    % spikeband = broadband2streamMinMax(filePrefix, blockId, outPrefix,options)

    %
    % outputs to a structure comparable to the continuous stream data structure
    % output format:
    %    
    %    spikeband.clock [1 x T]         :   clock times
    %    spikeband.minSpikeBand [96 x T] :   ms-binned min spikeband values for each channel
    %    spikeband.minSpikeInd [96 x T]  :   index of above
    %    spikeband.meanSquared           :   outputs MS value for a single channel, channel ID rotates ever 100 ms
    %    spikeband.meanSquaredChannel    :   channel index of above


if ~exist('options','var')
    options.foo = false;
end
if ~isstruct(options)
    error('broadband2streamMinMax: options must be a struct');
end

% filter?
if isfield(options,'filterType')
    filterType = options.filterType;
end
if ~exist('filterType','var')
    filterType = 'spikesmediumfiltfilt';
end

% use CAR? default is no
if ~isfield(options,'CAR')
    options.CAR = false;
end


% pick your spike Band filter
switch lower(filterType)
  case 'spikesmedium'
    filt = spikesMediumFilter();
  case 'spikeswide'
    filt = spikesWideFilter();
  case 'spikesnarrow'
    filt = spikesNarrowFilter();
  case 'spikesmediumfiltfilt'
    filt = spikesMediumFilter();
    useFiltfilt = true;
  case 'none'
    filt =[];
end
filt.PersistentMemory = false; % allow successive filtering


%% display information about the plx file
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information( inputfile );

disp(['Opened File Name: ' OpenedFileName]);
disp(['Version: ' num2str(Version)]);
disp(['Frequency : ' num2str(Freq)]);
disp(['Comment : ' Comment]);
disp(['Date/Time : ' DateTime]);
disp(['Duration : ' num2str(Duration)]);
disp(['Num Pts Per Wave : ' num2str(NPW)]);
disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);

% some information is only filled if the plx file version is >102
if ( Version > 102 )
    if ( Trodalness < 2 )
        disp('Data type : Single Electrode');
    elseif ( Trodalness == 2 )
        disp('Data type : Stereotrode');
    elseif ( Trodalness == 4 )
        disp('Data type : Tetrode');
    else
        disp('Data type : Unknown');
    end
        
    disp(['Spike Peak Voltage (mV) : ' num2str(SpikePeakV)]);
    disp(['Spike A/D Resolution (bits) : ' num2str(SpikeADResBits)]);
    disp(['Slow A/D Peak Voltage (mV) : ' num2str(SlowPeakV)]);
    disp(['Slow A/D Resolution (bits) : ' num2str(SlowADResBits)]);
end   

% there's some other info available. I don't really know what it means.
tmp = PL2GetFileIndex( inputfile );



%%% GLOBALS
% samplerate is kept in the Freq variable
NUM_SAMPLES_PLEXON_BLOCK = 65535; % PL2 data block size
MILLISECONDS_TO_GET = 120 * 1000;
% make that a multiple of the block size to simplify things
MILLISECONDS_TO_GET = ceil( MILLISECONDS_TO_GET / NUM_SAMPLES_PLEXON_BLOCK ) * NUM_SAMPLES_PLEXON_BLOCK;

% how many plexon samples per millisecond?
SAMPLES_PER_MS = Freq / 1000;
% make sure this is an int
assert( mod( SAMPLES_PER_MS, 1 ) == 0, 'not an integer number of samples' );

% Make data buffer this size
SAMPLES_IN_BUFFER = MILLISECONDS_TO_GET * SAMPLES_PER_MS;

% half of the "prior" buffer is to deal with edge effects
% the other half is to provide valid data for the previous frame's edge effects
EDGE_EFFECT_BUFFER = SAMPLES_PER_MS * 1000;
PRIOR_DATA_BUFFER_SAMPLES = EDGE_EFFECT_BUFFER * 2;
POST_DATA_BUFFER_SAMPLES = EDGE_EFFECT_BUFFER;

%% get some more info about the file
% 'slowcounts' tells us how many analog samples exist for each channel
[tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);
% find channels that have analog samples defined
analogChannels = find(slowcounts);
numChannels = numel( analogChannels );

% how many analog samples are there?
numSamples = unique( slowcounts( analogChannels ) );

% this better be the same for all channels, otherwise we don't know how to process.
assert( numel( numSamples ) == 1, 'different samples have different lengths. don''t know how to process' );

disp( sprintf('Num Analog Samples: %g', numSamples ) );

% keep some amount of previous data around, for edge effects
% create a buffer for prior data
priorData = zeros( PRIOR_DATA_BUFFER_SAMPLES, ...
                   numChannels, 'single');

% create a buffer to store moving avg of mean square value
MOV_AVG_BUFFER_SIZE = 100;
prevMsVals = zeros(MOV_AVG_BUFFER_SIZE, numChannels, 'single');



%% define the data span of interest
startIndex = 1;
endIndex = numSamples;
% make sure endIndex is an integer number of milliseconds
endIndex = floor( numSamples / SAMPLES_PER_MS ) * SAMPLES_PER_MS;

if isfield( options, 'startIndex' )
    startIndex = options.startIndex;
elseif isfield( options, 'startTime' )
    startIndex = floor( options.startTime * Freq ) + 1;
end
if isfield( options, 'endIndex' )
    endIndex = options.endIndex;    
elseif isfield( options, 'endTime' )
    endIndex = floor( options.endTime * Freq ) + 1;
end

%% may as well pre-allocate our millisecond-level representation
numMs = ( endIndex - startIndex + 1 ) / SAMPLES_PER_MS;
spikeband.clock = zeros( numMs, 1 );
spikeband.minSpikeBand = zeros( numMs, numChannels, 'single' );
spikeband.minSpikeBandInd = zeros( numMs, numChannels, 'uint8' );
spikeband.meanSquared = zeros( numMs, 1, 'single' );
spikeband.meanSquaredChannel = zeros( numMs, 1, 'uint16' );
spikeband.validSpikeBand = zeros( numMs*SAMPLES_PER_MS, numChannels, 'single' ); % FZ, for storing valid data
%spikeband.validSpikeBand = sparse( numMs*SAMPLES_PER_MS, numChannels ); % FZ, for dealing with OUT OF MEMORY error using zeros

% initialize some variables for our big for loop
currentBufferStartInd = 1;
numFramesProcessed = 0;
channelMsInd = 1;
lastBlockTime = 0;

% this is the buffer for filtering
% will be assembled from prior data and current data buffer
%dataForFiltering = zeros( numChannels, SAMPLES_IN_BUFFER + PRIOR_DATA_BUFFER_SAMPLES );
dataForFilter = zeros( SAMPLES_IN_BUFFER + PRIOR_DATA_BUFFER_SAMPLES, numChannels ); %fixed by FZ on 20191030

%% now process the data in an (awful, for now) loop
while currentBufferStartInd < endIndex
    tic
    disp( sprintf( 'Processing frame %02g, %02.1f %% done. Last frame took %02.1f seconds', ...
                   numFramesProcessed, ...
                   currentBufferStartInd / endIndex * 100, lastBlockTime ) );

    [ dataBuffer, samplesRead ] = readPL2Samples( inputfile, SAMPLES_IN_BUFFER, analogChannels, numFramesProcessed == 0 );
    %[ dataBuffer, samplesRead ] = readPL2Samples( inputfile, SAMPLES_IN_BUFFER, analogChannels - 128, numFramesProcessed == 0 ); % hard code to solve weird thing with wrong channel number

    % notch out noisy frequencies
    dataBuffer_normd = normalize(dataBuffer, 'centered');

    % set up a low pass filter to eliminate frequencies above 5000Hz % FZ added on 11/17/2019
    %Fs = 40000;
    %[b,a] = butter(4, [5000] / (Fs / 2), 'low');
    %dataBuffer_normd = filtfilt(b, a, dataBuffer_normd');
    %dataBuffer_normd = dataBuffer_normd';
    
    Fs = 40000;
    if sum(isnan(dataBuffer_normd(:))) == 0  % handle the last dataBuffer in the file, which has NaNs. Can't use pwelch to find peaks in this case.
        Pxx_normd = [];
        for nn = 1:size(dataBuffer_normd, 1)
            [ Pxx_normd(:,nn), w ] = pwelch( dataBuffer_normd( nn, : ), 40000*4 );
        end

        x = w * Fs/(2*pi);
        locs = {};
        pks = {};
        for nn = 1:size(dataBuffer_normd, 1)
            [pks_tmp, locs_tmp] = findpeaks(10 * log10(Pxx_normd(:, nn)'), x,  'MinPeakProminence', 5 );
            pksInBand = find(locs_tmp <= 5050 & locs_tmp >= 250); % FZ removed on 11/17/2019
            locs{nn} = locs_tmp(pksInBand); % FZ removed on 11/17/2019
                                            %locs{nn} = locs_tmp; % FZ added on 11/17/2019
            pks{nn} = pks_tmp(pksInBand); % FZ removed on 11/17/2019
                                          %pks{nn} = pks_tmp; % FZ added on 11/17/2019
        end
    end
    dataFiltered = dataBuffer_normd;
    parfor nn = 1:size(dataBuffer_normd, 1)

        % % dummy filter
        %hcas = dfilt.df1( 1, 1);
        for pk = 1:numel(locs{nn})
            freqToNotch = locs{nn}(pk); 
            %W0 = freqToNotch/(Fs/2); BW = 2/(Fs/2); %BW = W0/35;
            W0 = freqToNotch/(Fs/2); BW = 4/(Fs/2); %BW = W0/35; % FZ increased the bandwidth from 2Hz to 4Hz on 11/15/2019
            [b,a] = iirnotch(W0, BW);

            % cascade filters
            %hnew = dfilt.df1( b, a);
            %hcas = dfilt.cascade( hcas, hnew );

            % old filtering code
            dataFiltered(nn, :) = filter(b,a,dataFiltered(nn,:));
        end

        % % cascade filters
        %dataFiltered(nn, :) = hcas.filter( dataFiltered(nn,:) );
        
            %disp(['finished ' int2str(pk)])
    end


    %    keyboard
    dataBuffer = dataFiltered;

    % disp( sprintf( 'read %g / %g samples', samplesRead, SAMPLES_IN_BUFFER ) );

    currentBufferEndInd = currentBufferStartInd + samplesRead - 1;

    % check to make sure we have an integer number of milliseconds
    roundedNumSamples = SAMPLES_PER_MS * floor( samplesRead / SAMPLES_PER_MS );
    if roundedNumSamples ~= size( dataBuffer, 2)
        disp('broadband2StreamMin: dataBuffer unexpected size. OK if this is the end');
    end
    assert( mod( size( dataBuffer, 2) / SAMPLES_PER_MS, 1 ) == 0, 'read in # of samples is not integer multiple of 1 ms' );

    % convert to single (because faster? memory?)
    dataBuffer = single( dataBuffer' );

    if options.CAR
        % Apply common average referencing
        dataBuffer = dataBuffer - repmat( mean( dataBuffer, 2 ), 1, size( dataBuffer, 2) );
    end


    % concatenate prior data and most recent
    dataForFilter( 1 : PRIOR_DATA_BUFFER_SAMPLES, : ) = priorData;
    dataForFilter( PRIOR_DATA_BUFFER_SAMPLES + ( 1 : SAMPLES_IN_BUFFER ), : ) = dataBuffer;
    % store the end of most recent as new prior
    priorData = dataBuffer( end - PRIOR_DATA_BUFFER_SAMPLES + 1 : end, : );

    % now do the actual spikeband filtering!
    if ~isempty(filt)

        % FZ commented the follow if else statement on 190905
        %if useFiltfilt
        %    spikeBandData = filtfilthd( filt, dataForFilter );
        %else
        %    spikeBandData = filt.filter( dataForFilter );
        %end

        % FZ added the following on 190905
        filtLowCutoff = cutOffFreqs(1);
        filtHighCutoff = cutOffFreqs(2);
        [b,a] = butter(4, [filtLowCutoff filtHighCutoff] / (Fs / 2), 'bandpass');
        spikeBandData = single(filtfilt(b, a, double(dataForFilter)));
        % FZ modified above on 20191030 to solve memory issue
    else
        spikeBandData = single(dataForFilter);
    end

    % the beginning and end of this buffer are to be considered invalid
    validData = spikeBandData( EDGE_EFFECT_BUFFER + 1 : end-EDGE_EFFECT_BUFFER, : );

    %keyboard

    % there is also an EDGE_EFFECT_BUFFER length segment in the beginning
    % that was loaded with the last buffer, but has only been processed in
    % a valid way on this buffer
    validDataStartInd = currentBufferStartInd - EDGE_EFFECT_BUFFER;
    validDataEndInd = currentBufferEndInd - EDGE_EFFECT_BUFFER;
    % get the index of each ms-sample in terms of data sample indices
    sampleInds = validDataStartInd : SAMPLES_PER_MS : validDataEndInd;
    sampleInds_1 = validDataStartInd : 1 : validDataEndInd; % FZ for getting valid data
    keepInds = find( sampleInds > 0 ) : ...
        find( sampleInds <= endIndex, 1, 'last' );
    keepInds_1 = find( sampleInds_1 > 0 ) : ...
        find( sampleInds_1 <= endIndex, 1, 'last' ); % FZ for getting valid data

    % get the index of each ms-sample, but in ms
    sampleMsInds = ( (sampleInds-1) / SAMPLES_PER_MS ) + 1;

    indsToStore = sampleMsInds( keepInds );
    indsToStore_1 = sampleInds_1( keepInds_1 ); % FZ for getting valid data
    spikeband.clock( indsToStore ) = sampleInds( keepInds );
    
    % preallocate mins, max, moving avh
    mins = zeros( floor( size( validData )./ [ SAMPLES_PER_MS 1 ]),...
                  'single');
    minInd = zeros( floor( size( validData )./ [ SAMPLES_PER_MS 1 ]),...
                    'uint8');
    movAvgMs = zeros( floor( size( validData )./ [ SAMPLES_PER_MS 1 ]),...
                      'single');
    meanSquared = zeros( size( movAvgMs, 1), 1, 'single' );
    meanSquaredChannel = zeros( size( movAvgMs, 1), 1, 'uint8' );

    for nc = 1:size(validData,2)
        % get the spikeband min values and min indices
        cbroadband = reshape( validData( :, nc ), SAMPLES_PER_MS,[]);
        [channelMinVals channelMinInds] = min(cbroadband);
        % get min values per ms
        mins(:,nc) = channelMinVals;
        % and the sample index within the ms they occur
        minInd(:,nc) = uint8(channelMinInds);
        % get the mean-squared values
        localMs = mean(cbroadband.^2);
        
        latestMs = [prevMsVals(:,nc)' localMs];
        % take simple moving avg of prev ms points and current
        channelMovAvgMS=tsmovavg( latestMs, 's', MOV_AVG_BUFFER_SIZE );
        
        movAvgMs(:,nc) = channelMovAvgMS( MOV_AVG_BUFFER_SIZE + 1 : end );
        prevMsVals(:,nc) = localMs( end - MOV_AVG_BUFFER_SIZE + 1 : end );

    end

    % store the ms-level representations
    spikeband.minSpikeBand( indsToStore, : ) = mins( keepInds, : );
    spikeband.minSpikeBandInd( indsToStore, : ) = minInd( keepInds, : );
    spikeband.validSpikeBand( indsToStore_1, : ) = validData( keepInds_1, : );

    % store channel mean squared for only one channel at a time
    for i = 1:numel( keepInds )
        t = keepInds( i );
        channelMsInd = channelMsInd+1;
        if channelMsInd > numChannels, channelMsInd = 1; end
        meanSquared( t ) = movAvgMs( t, channelMsInd );
        meanSquaredChannel( t ) = uint8( channelMsInd );
    end
    spikeband.meanSquared( indsToStore ) = meanSquared( keepInds );
    spikeband.meanSquaredChannel( indsToStore ) = meanSquaredChannel( keepInds );
    
    
    % where should next buffer start
    currentBufferStartInd = currentBufferEndInd + 1;
    numFramesProcessed = numFramesProcessed + 1;
    lastBlockTime = toc;
end

%save
save(outputfile,'spikeband', '-v7.3');
%return 
clear all