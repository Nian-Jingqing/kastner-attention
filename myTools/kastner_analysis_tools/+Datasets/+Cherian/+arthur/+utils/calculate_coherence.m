function [ f, cxy ] = calculate_coherence( output, r_succ, new_sr, figs_dir )

% set n for fft options
n = 8;

% get original sample rate
orig_sr = r_succ.r( 1 ).sampleRate;

% get total number of trials
n_trials = numel( r_succ.r );
% get number of predictions
n_preds = numel( output );
% get number of channels
n_channels = size( output{ 1 }.full_true, 2 );

% get number of points for coherence
n_points = 2^(n+3)/2 + 1;

% intialize output tensor
cxy = zeros( n_preds, n_channels, n_points );

% for each prediction
for i_pred = 1:n_preds
    for i_channel = 1:n_channels
        % for each successful trial
        for i_trial = 1:n_trials
            % get original start and end indices
            orig_start_idx = r_succ.r( i_trial ).abstSidx;
            orig_end_idx = r_succ.r( i_trial ).abstEidx;
            
            % convert indices to new sample rate
            start_idx = round( orig_start_idx * ( new_sr / orig_sr ) );
            end_idx = round( orig_end_idx * ( new_sr / orig_sr ) );
            plot_range = start_idx:end_idx;

            % set x to true EMG and y to predicted EMG signal
            x = output{ i_pred }.full_true( plot_range, i_channel );
            y = output{ i_pred }.full_pred( plot_range, i_channel );

            % remove mean
            mean_rem_x = x - mean( x );
            mean_rem_y = y - mean( y );
            
            % calculate coherence between signals
            [ tmp_cxy, f ] = mscohere( mean_rem_x, mean_rem_y, 2^n, 2^(n-1), 2^(n+3), new_sr );
            % store in output tensor
            cxy( i_pred, i_channel, : ) = squeeze( cxy( i_pred, i_channel, : ) ) + tmp_cxy;
        end
    end
end

% average across all trials
cxy = cxy./n_trials;






