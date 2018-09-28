function [ in_data ] = resampleSpikes( mad, in_data, new_samplerate );
% get length and dim of input spikes
[ data_length, data_dim ] = size( in_data.spikes );

% get original sample rate
orig_samplerate = mad.analog(1).info.samplerate;

% get original number of points in analog data
orig_numpoints = mad.analog(1).info.numpoints;
resamp_numpoints = ceil( orig_numpoints * ( new_samplerate / orig_samplerate ) );
resamp_sampleperiod = 1 / new_samplerate;
resamp_bintimes = mad.analog( 1 ).info.starttime + ( 0 : resamp_numpoints - 1 )* resamp_sampleperiod;
resamp_binedges = resamp_bintimes - resamp_sampleperiod/2;
resamp_binedges( end + 1 ) = resamp_binedges( end ) + resamp_sampleperiod;
% find number of points at new sample rate and add extra bin to trim
%resamp_numbins = ceil( orig_numpoints * ( new_samplerate / orig_samplerate ) ) + 1;
% dont add extra bin since histcounts evenly splits 

%numBins = resamp_numbins;
for nc = 1:data_dim
    % get binned spikes
    orig_spike_inds = find( in_data.spikes( :, nc ) == 1 );
    % get orig spike times
    orig_spike_times = orig_spike_inds / orig_samplerate;
    binnedSpks = histc( orig_spike_times, resamp_binedges );
    % push spikes to output struct
    in_data.rebinned_spikes( :, nc ) = binnedSpks( 1 : end-1 );
end
