function [ madSpikes, spikeInfo ] = extractSpikes( mad, excludeN, sampleRate )
% extract spiking data from madStruct

% mad : madStruct
% excludeN : neurons to exclude from output spike matrix
% sampleRate : sample rate of analog data

% madSpikes : column matrix of spiking data ( time x nchannels )
% spikeInfo : row matrix containing the information about cluster/channel
    
% number of neural channels
    nChannels = size( mad.spikes, 1 );
    % number of neural channels to remove from output matrix
    numRemoveNeuralChan = size( excludeN, 2 );
    % number of samples in output matrix
    outputDataSize = numel( mad.analog( 1 ).waveform );
    % analog offset time
    analogStartTime = mad.analog( 1 ).info.starttime;
    % intialize output spike matrix
    madSpikes = zeros( nChannels, outputDataSize );    
    % initialize tmp matrices with same size as output matrix
    % w/ 3 extra (channel #/cluster #/original channel #)    
    tmpSpikes = zeros( nChannels, outputDataSize+3 );
    madSpikesSort = zeros( nChannels, outputDataSize+3 );
    
    % extract Spikes from MadStruct and push to output spike matrix
    
    for i = 1 : size( mad.spikes, 1 )
        % extract original timestamps
        tstamps_init = mad.spikes( i ).timestamps;
        % offset spiking timestamps to be aligned with the analog channel timing
        tstamps_offset = tstamps_init - analogStartTime;
        % convert timestamps to sample index numbers that align with analog data
        tstamps = Rbuild.Utils.convert2sampleNumber( ...
            tstamps_offset, sampleRate );
        % remove spikes that happen before analog channel recording started        
        spike_idx = tstamps( tstamps > 0 & tstamps <= outputDataSize );
        % transpose to row vector
        spike_idx = spike_idx';
        % initialize vector for channel binned spiking data
        spikes = zeros( 1, outputDataSize );
        % set spike events to 1
        spikes( spike_idx ) = 1;
        % extract info about spike and channel and cluster
        channel = mad.spikes( i ).info.channel;
        cluster = mad.spikes( i ).info.cluster;
        % push information to the output matrix
        tmpSpikes( i, 1 ) = channel;
        tmpSpikes( i, 2 ) = cluster;
        tmpSpikes( i, 3 ) = i;
        tmpSpikes( i, 4:end ) = spikes;
    end
    
    % organize spikes matrix in ascending order of channel/cluster number
    madSpikesSort = sortrows( tmpSpikes );
    % remove specified excluded spiking channels
    if ~isempty( excludeN )
            madSpikesSort( excludeN, : ) = [];
    end
    assert( size( madSpikes, 2 ) == ...
            size( madSpikesSort( :, 4:end ), 2 ), ...
            'ERROR: madSpikesSort does not match dimensionality of output analog data.');    
    madSpikes = madSpikesSort( :, 4:end )';
    spikeInfo = madSpikesSort( :, 1:3 );
end
