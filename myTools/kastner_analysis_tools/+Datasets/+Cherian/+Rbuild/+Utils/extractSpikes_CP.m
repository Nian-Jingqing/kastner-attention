function [ data ] = extractSpikes_CP( mad )

% get info about analog channels
analog_samplerate = mad.analog(1).info.samplerate;
analog_sampleperiod = 1/ analog_samplerate;
analog_bintimes = mad.analog( 1 ).info.starttime + ( 0 : mad.analog( 1 ).info.numpoints - 1 )* analog_sampleperiod;

analog_binedges = analog_bintimes - analog_sampleperiod/2;
analog_binedges( end + 1) = analog_binedges( end ) + analog_sampleperiod;

%% get the spike times
binEdges = analog_binedges;
for nc = 1:numel( mad.spikes )
    binnedSpks = histc( mad.spikes( nc ).timestamps, binEdges );
    % last item is reserved for times exactly matching bin edge, should be able to trim
    data.spikes( :, nc ) = binnedSpks( 1: end-1 );
end