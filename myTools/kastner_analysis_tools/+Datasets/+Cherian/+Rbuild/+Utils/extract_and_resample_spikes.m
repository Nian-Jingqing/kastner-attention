function [ data ] = extract_and_resample_spikes( mad, new_samplerate )


try
    % get original sample rate
    analog_samplerate = mad.analog( 1 ).info.samplerate;
    % get original number of points in analog data
    analog_numpoints = mad.analog( 1 ).info.numpoints;
    % get analog start time
    analog_starttime = mad.analog( 1 ).info.starttime;
catch ME
    analog_samplerate = mad.analoginfo( 1 ).samplerate;
    analog_numpoints = mad.analoginfo( 1 ).numpoints;
    analog_starttime = mad.analoginfo( 1 ).starttime;
end

% calculate number of points at new sample rate
resamp_numpoints = ceil( analog_numpoints * ( new_samplerate / analog_samplerate ) );
% calculate sample period
resamp_sampleperiod = 1 / new_samplerate;
% create vector of bin times
resamp_bintimes = analog_starttime + ( 0 : resamp_numpoints - 1 )* resamp_sampleperiod;
% offset by half the period to get the bin edges
resamp_binedges = resamp_bintimes - resamp_sampleperiod/2;
% add final bin at end for right bin edge (will be trimmed from output)
resamp_binedges( end + 1 ) = resamp_binedges( end ) + resamp_sampleperiod;

% for each neuron
% intialize output matrix
%spikes = zeros( resamp_numpoints, numel( mad.spikes ) );
for nc = 1:numel( mad.spikes )
    
    % extract info about spike and channel and cluster
    try
        channel = mad.spikes( nc ).info.channel;
        cluster = mad.spikes( nc ).info.cluster;
    catch ME
        channel = mad.spikesinfo( nc ).channel;
        cluster = mad.spikesinfo( nc ).cluster;
    end
    
    % get binned spikes
    binned_spks = histc( mad.spikes( nc ).timestamps, resamp_binedges );
    % push spikes to output struct
    data.spikes( :, nc ) = binned_spks( 1 : end-1 );
    data.spikeInfo( nc, : ) = [ channel cluster ];
end
