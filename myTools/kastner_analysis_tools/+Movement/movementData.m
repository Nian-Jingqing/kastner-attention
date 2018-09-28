classdef movementData < R.Rstruct   
    methods
        % constructor
        function obj = movementData( varargin )
            obj = obj@R.Rstruct (varargin{:} );
        end
        function [ ] = find_align_point_times( r, speed_threshold_perc, align_field_name )
        % Written by LW
            % find total number of trials in Rstruct
            n_trials = size( r.r, 2 );
            for i = 1:n_trials
                % calculate speed
                s = sqrt( sum( r.r(i).kin(3:4, :).^2 ) );
                
                % find index and value for max speed
                [max_s, max_s_idx] = max( s );

                % calculate threshold for movement onset based on percentage of max speed
                s_threshold = max_s * speed_threshold_perc;

                % find alignment point index
                align_idx = find( s >= s_threshold, 1 );

                % store value in rstruct
                r.r( i ).( align_field_name ) = align_idx;
             end
        end
        function [ keep_trials ] = find_align_chop_keep_trials( r, ap_threshold, max_speed_min, pre_post_ap_ms, align_field_name )
        % Written by LW
            pre_ap_ms = pre_post_ap_ms( 1 );
            post_ap_ms = pre_post_ap_ms( 2 );
            keep_trials = false( 1, numel( r.r ) );
            n_trials = numel( keep_trials );
            %align_field_name = 'move_align';
            r.find_align_point_times( ap_threshold, align_field_name );
            
            for i= 1:n_trials
                % if trial is success
                if ( r.r( i ).trialType + 0 ) == 1
                    % calculate speed for trial
                    s = sqrt( sum( r.r(i).kin(3:4, :).^2 ) );
                    % find value and index for max speed
                    [ max_s, max_s_idx ] = max( s );
                    % plot trial if trace meets conditions
                    if max_s >= max_speed_min
                        % get movement onset idx
                        align_idx = r.r( i ).( align_field_name );
                        try
                            % if test fails then trial is not kept
                            test = s( align_idx - pre_ap_ms : align_idx + post_ap_ms - 1 );
                            keep_trials( i ) = true;
                        end % try
                    end % if max_s
                end % if r.r
            end % for
        end % find_keep_trials
        function [ keep_trials ] = find_overlap_chop_keep_trials( r )
        % Written by LW
        % currently just returns all trials but made for future potential implementation
            keep_trials = 1:numel( r.r );
        end % find_olap_chop_keep_trials
        
        function [ data ] = generate_align_chop_lfads_data( r, pre_post_ap_times, keep_trials, exclude_channels, align_field_name )
        % Written by LW
        % assumes 1 ms binned data in R struct
            pre_ap_ms = pre_post_ap_times( 1 );
            post_ap_ms = pre_post_ap_times( 2 );
            
            % trial length
            n_samples = numel( -1*pre_ap_ms:post_ap_ms ) - 1;
            
            % number of tot trials
            n_tot_trials = size( keep_trials, 2 );

            % number of keep trials
            n_keep_trials = sum( keep_trials );

            % find indices of keep trials
            keep_trial_inds = find( keep_trials == 1 );
            
            % original number of channels
            tot_channels = size( r.r(1).spikes, 1 );
            % number of channels after high xcorr removal
            n_channels = tot_channels - numel( exclude_channels );

            % intialize logical vector to index spikes
            keep_channels = true( 1, tot_channels );
            % set exclude channels for removal
            keep_channels( exclude_channels ) = false;
            
            % initialize output counts matrix
            counts = zeros( n_keep_trials, n_channels, n_samples );
            % intialize output condition Id vector
            condition_id = zeros( 1, n_keep_trials );
            % intialize output timeVecMs
            tmp = 1:( numel( -1*pre_ap_ms:post_ap_ms ) - 1 );
            time_vec_ms = tmp - pre_ap_ms;
            
            for i=1:n_keep_trials
                i_trial = keep_trial_inds( i );
                % get spikes
                spikes = r.r( i_trial ).spikes;
                trial_cond = r.r( i_trial ).conditionId;
                % get alignment point idx
                align_idx = r.r( i_trial ).( align_field_name );
                % get start and end idx for trials
                start_idx = align_idx - pre_ap_ms;
                end_idx = align_idx + post_ap_ms - 1;
                counts( i, :, : ) = spikes( keep_channels, start_idx:end_idx );
                condition_id( i ) = trial_cond;
            end
            data = struct();
            data.counts = counts;
            data.timeVecMs = time_vec_ms;
            data.conditionId = condition_id;
        end
        function [ data ] = generate_overlap_chop_lfads_data( r, keep_trials, trial_time_ms, trial_olap_ms, exclude_channels )
            
            % get spikes in continous form
            spikes = [ r.r( keep_trials ).spikes ];
            % remove excluded channels
            spikes( exclude_channels, : ) = [];
            % get length of data in ms (pretty easy since 1ms binned data)
            data_time_ms = size( spikes, 2 );
            % determine number of channels
            n_dim = size( spikes, 1 );
            % determine number of trials
            n_trials = ceil( ( data_time_ms - trial_time_ms ) / ( trial_time_ms - trial_olap_ms ) ) + 1;

            % intialize counts
            counts = zeros( n_trials, n_dim, trial_time_ms );
            % now generate counts
            for i_trial=1:n_trials
                % get start idx (ms)
                start_idx = ( i_trial-1 )*trial_time_ms - ( i_trial-1 )*trial_olap_ms + 1;
                % get end idx( ms)
                end_idx = ( i_trial-1 )*trial_time_ms - ( i_trial-1 )*trial_olap_ms + trial_time_ms;
                % try to extract trial
                try
                    counts( i_trial, :, : ) = spikes( :, start_idx:end_idx );
                % but if you go over index limit of matrix
                catch ME
                    % start from the end of the data and have more overlap for the last trial
                    start_idx = data_time_ms - trial_time_ms + 1;
                    end_idx = data_time_ms;
                    counts( i_trial, :, : ) = spikes( :, start_idx:end_idx );
                end
            end
            % intialize output struct
            data = struct();
            data.counts = counts;
            data.timeVecMs = 1:trial_time_ms;
            data.conditionId = 1:n_trials;
        end
        
        function[ factors, rates, gen_states, con_inputs ] = get_output_from_lfads( r, run, ds )
            disp( 'INFO: Getting spiking data from datasets ...' )
            for i=1:numel( run.params.i_block )
                block = run.params.i_block( i );
                data = ds.loadData( block, run.params.dataset_type );
                if i == 1
                    spikes = [ data.r.spikes ];
                else
                    spikes = [ spikes [ data.r.spikes ] ];
                end
            end
            
            % load output of run
            disp( 'INFO: Loading sequence data and posterior means ...' )
            run.loadSequenceData();
            run.loadPosteriorMeans();
            run.addPosteriorMeansToSeq();

            % get factors and rates form posterior means
            tot_factors = run.posteriorMeans.factors;
            tot_rates = run.posteriorMeans.rates;
            % get generator states
            tot_gen_states = [ run.sequenceData{ 1 }.generator_states ];
            % get controller inputs
            tot_con_inputs = [ run.sequenceData{ 1 }.controller_outputs ];
            % get geneator dim and individual trial length
            [ gen_dim, trial_length ] = size( run.sequenceData{ 1 }(1).generator_states );
            % get controller dim
            con_dim = size( tot_con_inputs, 1 );

            % determine number of trials
            n_trials = numel( run.sequenceData{ 1 } );
            tot_gen_states = reshape( tot_gen_states, [ n_trials, gen_dim, trial_length ] );
            tot_con_inputs = reshape( tot_con_inputs, [ n_trials, con_dim, trial_length ] );

            % permute data to trials x channels x time
            tot_factors = permute( tot_factors, [ 3 1 2 ] );
            tot_rates = permute( tot_rates, [ 3 1 2 ] );

            trialize_method = run.params.trialize_method;
            
            if strcmp( trialize_method, 'align-chop' ) % if trials generated using 'align-chop'
                factors = tot_factors;
                rates = tot_rates;
                
            else
                            
                disp( 'INFO: Extracting orignal data from overlapped trials ...' )
                % Extract relevant hyperparameters
                trialTimeMs = run.params.trial_time_ms;
                trialOlapMs= run.params.trial_olap_ms;
                binWidthMs = run.params.spikeBinMs;
                
                load( ds.exclude_neuron_path )
                % transpose spike to column matrix
                spikes = spikes';
                
                spikes( :, all_neurons_to_remove ) = [];
                
                nbins = ceil( size( spikes, 1 ) / binWidthMs );
                
                % Find trial length parameters in terms of data sampling frequency
                trialTime = trialTimeMs / binWidthMs;
                trialOlap = trialOlapMs / binWidthMs;
                
                % Find dimensionality of data
                nFactors = size( tot_factors, 2 );
                nSpikes = size( spikes, 2 );
                nRates = size( tot_rates, 2 );
                nGen = size( tot_gen_states, 2 );
                nCon = size( tot_con_inputs, 2 );
                %nbins = size( spikes, 1 );
                
                % Dimensionality Check
                %assert(nSpikes == nRates,...
                %       'ERROR: Dimensionality of spiking data and LFADS inferred rates do not match.')
                
                % Intialize output data matrices
                factors = zeros( nFactors, nbins );
                rates = zeros( nRates, nbins );
                gen_states = zeros( nGen, nbins );
                con_inputs = zeros( nCon, nbins );
                % Design filters
                lfadsFilter1 = r.designTrialFilter( nbins, trialTime, trialOlap, 0 );
                lfadsFilter2 = r.designTrialFilter( nbins, trialTime, trialOlap, 1 );
                lfadsFilter3 = r.designTrialFilter( nbins, trialTime, trialOlap, 2 );
                
                % Extract number of trials 
                nTrials = size( tot_factors, 1 );

                for i=1:nTrials
                    % Choose filter
                    %disp(i)
                    if i == 1
                    lfadsTrialFilter = lfadsFilter1;
                    elseif i == size(tot_factors, 1 )
                        lfadsTrialFilter = lfadsFilter3;
                    else
                        lfadsTrialFilter = lfadsFilter2;
                    end
                    
                    % Extract trial data
                    rawTrialFactors = squeeze( tot_factors( i, :, : ) );
                    rawTrialRates = squeeze( tot_rates( i, :, : ) );
                    rawTrialGen = squeeze( tot_gen_states( i, :, : ) );
                    rawTrialCon = squeeze( tot_con_inputs( i, :, : ) );         
                    % Filter trial data
                    filtTrialFactors = r.filterChannel( rawTrialFactors, lfadsTrialFilter);
                    filtTrialRates = r.filterChannel( rawTrialRates, lfadsTrialFilter );
                    filtTrialGen = r.filterChannel( rawTrialGen, lfadsTrialFilter );
                    filtTrialCon = r.filterChannel( rawTrialCon, lfadsTrialFilter );
                    
                    if i~= size( tot_factors, 1 )
                        trialStart = 1 + (i-1)*trialTime - (i-1)*trialOlap;
                        trialEnd = trialStart + trialTime - 1;
                    else
                        trialStart = nbins - trialTime + 1;
                        trialEnd = nbins;
                    end
                    
                % Add trial data to output channels
                factors( :, trialStart:trialEnd ) = factors( :, trialStart:trialEnd ) + filtTrialFactors;
                rates( :, trialStart:trialEnd ) = rates( :, trialStart:trialEnd ) + filtTrialRates;
                gen_states( :, trialStart:trialEnd ) = gen_states( :, trialStart:trialEnd ) + filtTrialGen;
                con_inputs( :, trialStart:trialEnd ) = con_inputs( :, trialStart:trialEnd ) + filtTrialCon;
                
                end
                factors = factors';
                rates = rates';
                gen_states = gen_states';
                con_inputs = con_inputs';
            end
        end
        function [ filt_main ] = designTrialFilter( r, nbins, trialTime, trialOlap, filterType )
        % Design ascending linear filter
            leftLinFilt = 1:trialOlap;
            leftLinFilt = leftLinFilt/trialOlap;
            % Design descending linear filter
            rightLinFilt = trialOlap-1:-1:0;
            rightLinFilt = rightLinFilt/trialOlap;
            
            switch filterType
              case 0
                % Filter for first trial
                filt_main = ones(1,trialTime);
                filt_main(end-trialOlap+1:end) = rightLinFilt;
              case 1
                % Filter for middle trials                
                filt_main = ones(1,trialTime);
                filt_main(1:trialOlap) = leftLinFilt;
                filt_main(end-trialOlap+1:end) = rightLinFilt;
              case 2
                % Filter for last trial
                trialNewData = mod(nbins-trialTime,trialTime-trialOlap);
                nTrialIgnore = trialTime - trialNewData - trialOlap;
                trialIgnore = 1:nTrialIgnore;
                trialFilter = nTrialIgnore:nTrialIgnore+trialOlap-1;
                filt_main = ones(1,trialTime);
                filt_main(trialIgnore) = 0;
                filt_main(trialFilter) = leftLinFilt;
            end
        end
        function [ filtChannel ] = filterChannel( r, channel, filter )
            filtChannel = channel .* filter;
        end
        function [ trial_start_inds, trial_end_inds ] = get_data_inds_lists( r, data_switch )
        % get trial start and trial end indices in two separate vectors
        % this function's output will get passed to `gen_data_inds`
        
        % get trial type enumerator and typecast to integer by adding 0
            trial_type = [ r.r.trialType ] + 0;
            

            % create logical vector for trials to keep in Rstruct
            keep_trials = false( 1, numel( r.r ) ); 

            switch data_switch
              case 0 % all data (trial/intertrial)
                keep_trials( : )  = true;
              case 1 % successful trials only
                keep_trials(find( trial_type == 1 )) = true;
              case 2 % all trials
                keep_trials(find( trial_type ~= 0 )) = true;
              case 3 % all intertrial
                keep_trials(find( trial_type == 0 )) = true;
            end
            % extract keep trials
            r_keep = r.copy();
            r_keep.r = r_keep.r( keep_trials );

            % get trial indices
            trial_start_inds = [ r_keep.r.abstSidx ];
            trial_end_inds = [ r_keep.r.abstEidx ];
        end
        function [ data_inds ] = gen_data_inds( r, run, data_switch, new_data_length, new_sr )
        % generate logical vector to extract relevant data from raw_spikes
        % intialize vector to hold all sample ranges

            par = run.params;
            orig_sr = r.r(1).sampleRate;
            %new_sr = 1000 / par.spikeBinMs;
            
            n_trial_samples = round( ( par.trial_time_ms/1000 ) * new_sr );
            
            % get data indices at orig_sr
            [ orig_start_inds, orig_end_inds ] = r.get_data_inds_lists( data_switch );

            % resample down to new_sr
            new_start_inds = round( orig_start_inds * ( new_sr / orig_sr ));
            new_end_inds = round( orig_end_inds * ( new_sr / orig_sr ));
            
            % fix indices less than 1 or greater than length of new_data_length
            new_start_inds( new_start_inds < 1 ) = 1;
            new_end_inds( new_end_inds > new_data_length ) = new_data_length;

            % make sure same number of start and end inds
            assert( numel( new_start_inds ) == numel( new_end_inds ), ...
                    'ERROR: trial_start and trial_end dont match. ' )
            tmp_data_inds = [];
            % add all success ranges to single vector
            for i=1:numel( new_start_inds )
                tmp_vec_inds = new_start_inds( i ):new_end_inds( i );
                tmp_data_inds = [ tmp_data_inds tmp_vec_inds ];
            end
            % intialize output logical vector
            data_inds = false( 1, new_data_length );
            data_inds( tmp_data_inds ) = true;
        end
        
        function [ ] = plot_avg_emg( r, pre_post_mo_ms, plot_title, channel )
        % Written by LW
        % this function isnt very helpful right now. Need to have some sort of phase alignment step fo the emg signals to get meaningful average
        % initialize color map
            cmap = hsv;
            n_conds = numel( unique( [ r.r.conditionId ] ) );
            
            % determine number of trials
            n_trials = size( r.r, 2 );

            % extract pre movement onset and post movement onset times
            pre_mo_ms = pre_post_mo_ms( 1 );
            post_mo_ms = pre_post_mo_ms( 2 );

            tmp = 1:( numel( -1*pre_mo_ms:post_mo_ms ) - 1 );
            tmp = tmp - pre_mo_ms;
            
            % intialize figure
            figure( 404 )

            for i_cond = 1:n_conds
                % get color for trace
                clrind = floor( (i_cond - 1) / n_conds * size(cmap, 1) ) + 1;
                
                % get indices of trials for condition
                cond_trials = find( [ r.r.conditionId ] == i_cond );
                n_samples = numel( -1*pre_mo_ms:post_mo_ms ) - 1;
                sum_emg = zeros( 1, n_samples );
                for i_trial = 1:numel( cond_trials )
                    move_onset_idx = r.r( i_trial ).moveOnset;
                    start_idx = move_onset_idx - pre_mo_ms;
                    end_idx = move_onset_idx + post_mo_ms - 1;
                    
                    % get trial emg
                    trial_emg = r.r( i_trial ).lp_emg( channel, start_idx: end_idx );
                    if i_trial > 1
                        prev_trial_emg = r.r( i_trial-1 ).lp_emg( channel, start_idx: end_idx );
                        delay = finddelay( prev_trial_emg, trial_emg );
                        trial_emg = r.r( i_trial ).lp_emg( channel, start_idx+delay: end_idx+delay );
                    end
                    % add 
                    sum_emg = sum_emg + trial_emg;
                end
                % average the summed emg signal for the condition across the number of trials for cond
                avg_emg = sum_emg / numel( cond_trials );
                h = plot( tmp, avg_emg );
                set( h, 'Color', cmap( clrind, : ) );
                set( h, 'LineWidth', 4 );
                hold on
                %if i_cond == 1
                %    break
                %end
            end
            axis( 'tight' )
            title( plot_title )
            xlabel( 'Time Aligned to Move Onset (ms)' )
        end
        function [ resamp_data ] = resample_data( r, data, orig_sr, new_sr )
            % if data is organized as a row matrix, transpose
            if size( data, 1 ) < size( data, 2 )
                data = data';
            end
            
            resamp_data = resample( data, new_sr, orig_sr );
        end
        function [ rebin_data ] = rebin_data( r, data, new_sr )
        % get data
            bin_size = 1000 / new_sr;
            % if data is organized as a column matrix, transpose
            if size( data, 1 ) > size( data, 2 )
                data = data';
            end
            % get data dim
            data_dim = size( data, 1 );
            % get data length
            data_length = size( data, 2 );
            % find number of bins
            n_bins = floor( data_length / bin_size );

            points_to_keep = n_bins * bin_size;
            data = data( :, 1:points_to_keep );

            data2 = reshape( full( data ), data_dim, bin_size, n_bins );
            rebin_data = squeeze( sum( data2, 2 ) );

            if data_dim == 1
                rebin_data = rebin_data( : )';
            end
        end
        function [ resamp_data ] = resample_field( r, run, data_field )
        % get data
            data = [ r.r.( data_field ) ];
            % get original sample rate
            orig_sr = r.r(1).sampleRate;
            % get new sample rate
            new_sr = 1000 / run.params.spikeBinMs;

            % if data is organized as a row matrix, transpose
            if size( data, 1 ) < size( data, 2 )
                data = data';
            end
            
            resamp_data = resample( data, new_sr, orig_sr );
        end
        
        function [ rebin_data ] = rebin_field( r, run, data_field )
        % get data
            data = [ r.r.( data_field ) ];
            % get bin size
            bin_size = run.params.spikeBinMs;
            
            % if data is organized as a column matrix, transpose
            if size( data, 1 ) > size( data, 2 )
                data = data';
            end
            % get data dim
            data_dim = size( data, 1 );
            % get data length
            data_length = size( data, 2 );
            % find number of bins
            n_bins = ( data_length / bin_size );

            points_to_keep = n_bins * bin_size;
            data = data( :, 1:points_to_keep );

            data2 = reshape( full( data ), data_dim, bin_size, n_bins );
            rebin_data = squeeze( sum( data2, 2 ) );

            if data_dim == 1
                rebin_data = rebin_data( : )';
            end
        end
        function [ ] = old_rebin_field( r, infield, outfield, binSize )
        % iterate over all trials
            for itrial = 1:numel( r.r )

                % iterate over all the fields we want to bin
                d = r.r( itrial ).( infield );
                dataDim = size( d, 1);
                
                % trim any data that won't fit in an even number of bins
                dataLength = size( d, 2);
                numBins = floor( dataLength / binSize );
                pointsToKeep = numBins * binSize;
                d = d( :, 1:pointsToKeep );
                
                % reshape, sum, and squeeze
                d2 = reshape( full( d ), dataDim, binSize, numBins );
                dataBinned = squeeze( sum( d2, 2 ) );
                
                % if the original data is 1-dimensional, there are some problems with squeeze
                %  correction for that is below
                if dataDim == 1
                    dataBinned = dataBinned(:)';
                end
                
                % assign that data
                r.r( itrial ).( outfield ) = dataBinned;
            end % itrial
        end % binData

        function [ ] = plot_trial_emg( r, pre_post_mo_ms, plot_title, channel, plot_conds )
        % Written by LW
        % initialize color map
            cmap = hsv;
            numConds = numel( unique( [ r.r.conditionId ] ) );
            assert( max( plot_conds ) <= numConds, 'ERROR: Plotting nonexistent conditions. [ Invalid range ]' )
            % determine number of trials
            n_trials = size( r.r, 2 );

            % extract pre movement onset and post movement onset times
            pre_mo_ms = pre_post_mo_ms( 1 );
            post_mo_ms = pre_post_mo_ms( 2 );

            % intialize figure
            f = Plot.blankFigure( 404 );
            clf;
            
            opacity = 0.8;
            lineStyle = '-';
            lineWidth = 3.5;

            % for all trials
            for i_trial=1:n_trials
                % get condition id
                icond = r.r( i_trial ).conditionId;
                if ~isempty( find( plot_conds == icond ) )
                    % get color for trace
                    clrind = floor( (icond - 1) / numConds * size(cmap, 1) ) + 1;
                    % get move onset
                    move_onset_idx = r.r( i_trial ).moveOnset;
                    % find start index
                    start_idx = move_onset_idx - pre_mo_ms;
                    % find end index
                    end_idx = move_onset_idx + post_mo_ms - 1;
                    % temp plot range
                    tmp = 1:( numel( -1*pre_mo_ms:post_mo_ms ) - 1 );
                    % adjusted plot range for move onset
                    plot_range = tmp - pre_mo_ms;
                    % get trial emg for channel
                    trial_emg = r.r( i_trial ).lp_emg( channel, start_idx: end_idx );
                    %h = plot( tmp, trial_emg );
                    % plot emg
                    h = Plot.patchline( plot_range , trial_emg, ...
                                        'LineWidth', lineWidth, ...
                                        'LineStyle', lineStyle, ...
                                        'edgeAlpha', opacity, ...
                                        'edgeColor', cmap( clrind, : ), ...
                                        'faceAlpha', 0.0, ...
                                        'faceColor', cmap( clrind, : ), ...
                                        'AlignVertexCenters', 'off' );
                    hold on
                end
            end
            axis( 'tight' )
            title( plot_title )
            xlabel( 'Time Aligned to Move Onset (ms)' )
            ylabel( 'EMG' )
        end
        function [ trial_conds, keep_trials ] = plot_speed( r, pre_post_mo_ms, plot_title, max_speed_min )
        % Written by LW
        % extract pre movement onset and post movement onset times
            pre_mo_ms = pre_post_mo_ms( 1 );
            post_mo_ms = pre_post_mo_ms( 2 );

            % determine number of trials
            n_trials = size( r.r, 2 );
            
            % intialize figure
            figure( 303 );
            clf;
            % intialize count
            count = 1;
            
            % set x vector to be aligned to movement onset
            tmp = 1:( numel( -1*pre_mo_ms:post_mo_ms ) - 1 );
            plot_range = tmp - pre_mo_ms;

            % initialize color map
            cmap = hsv;
            numConds = numel( unique( [ r.r.conditionId ] ) );
            keep_trials = false( 1, numel( r.r ) );

            % plot settings
            opacity = 0.6;
            lineStyle = '-';
            lineWidth = 2.5;

            % for each trial
            for i = 1: n_trials
                % calculate speed
                s = sqrt( sum( r.r(i).kin(3:4, :).^2 ) );
                % get condition id
                icond = r.r( i ).conditionId;
                % get color for trace
                clrind = floor( (icond - 1) / numConds * size(cmap, 1) ) + 1;
                % get move onset
                move_onset_idx = r.r( i ).moveOnset;
                % find value and index max speed
                [max_s, max_s_idx] = max( s );
                try
                    if max_s >= max_speed_min 
                        % plot speed trace
                        %h = plot( tmp, s( move_onset_idx - pre_mo_ms: move_onset_idx + post_mo_ms - 1) );
                        h = Plot.patchline( plot_range , s( move_onset_idx - pre_mo_ms: move_onset_idx + post_mo_ms - 1), ...
                                        'LineWidth', lineWidth, ...
                                        'LineStyle', lineStyle, ...
                                        'edgeAlpha', opacity, ...
                                        'edgeColor', cmap( clrind, : ), ...
                                        'faceAlpha', 0.0, ...
                                        'faceColor', cmap( clrind, : ), ...
                                        'AlignVertexCenters', 'off' );
                        %set( h, 'Color', cmap( clrind, : ) );
                        hold on
                        % add condition to keep trial condition vector
                        trial_conds( count ) = icond;
                        % change keep trial flag to true
                        keep_trials( i ) = true;
                        count = count + 1;
                    else
                        disp( 'trial max speed too low' )
                    end
                catch ME
                    disp( ME.identifier )
                    disp( 'trial too short.' )
                end
            end
            axis( 'tight' )
            ylabel( 'Speed (cm/s)' )
            xlabel( 'Time Aligned to Move Onset (ms)' )
            title( plot_title )
        end

        function [ ] = find_move_onset_times( r, speed_threshold_perc )
        % Written by LW
        % LW - removed argument go_cue_field from function since it wasnt used
        % find total number of trials in Rstruct
            n_trials = size( r.r, 2 );
            for i = 1:n_trials
                % extract x and y velocities
                %v_x = r.r( i ).kin( 3, : );
                %v_y = r.r( i ).kin( 4, : );
                
                % calculate speed
                %s = sqrt( v_x.^2 + v_y.^2 );
                s = sqrt( sum( r.r(i).kin(3:4, :).^2 ) );
                
                % find max speed
                % find index for max speed, move onset, and go cue
                [max_s, max_s_idx] = max( s );


                %go_cue = r.r( i ).( go_cue_field );
                %s( 1:go_cue ) = 0;

                % find max speed during trial
                %max_speed = max( s );
                % calculate threshold for movement onset based on percentage of max speed
                s_threshold = max_s * speed_threshold_perc;

                % find move onset and store value in Rstruct
                move_onset_idx = find( s >= s_threshold, 1 );
                r.r( i ).moveOnset = move_onset_idx;

                tmp = 1:size(r.r(i).kin,2);
                tmp = tmp - r.r( i ).relGoCue;

                %plot(tmp, s);
                %hold on;
                %plot(move_onset_idx - r.r( i ).relGoCue, s(move_onset_idx),'o');

                %                plot( s );
                %                hold on
                %if move_onset_idx - r.r( i ).relGoCue < 50
                    %    keyboard
                    %end
            end
        end
        
        function [ keep_trials ] = find_keep_trials( r, mo_threshold, max_speed_min, pre_post_mo_ms )
        % Written by LW
            pre_mo_ms = pre_post_mo_ms( 1 );
            post_mo_ms = pre_post_mo_ms( 2 );
            keep_trials = false( 1, numel( r.r ) );
            n_trials = numel( keep_trials );
            r.find_move_onset_times( mo_threshold );
            for i= 1:n_trials
                % if trial is success
                if ( r.r( i ).trialType + 0 ) == 1
                    % calculate speed for trial
                    s = sqrt( sum( r.r(i).kin(3:4, :).^2 ) );
                    
                    % find value and index for max speed
                    [ max_s, max_s_idx ] = max( s );

                    % plot trial if trace meets conditions
                    if max_s >= max_speed_min
                        % get movement onset idx
                        move_on_idx = r.r( i ).moveOnset;
                        try
                            % if test fails then trial is not kept
                            test = s( move_on_idx - pre_mo_ms : move_on_idx + post_mo_ms - 1 );
                            keep_trials( i ) = true;
                        end % try
                    end % if max_s
                end % if r.r
            end % for
        end % find_keep_trials

        function [ data ] = generate_lfads_data_2( r, pre_post_mo_times, keep_trials, fieldToGen, alignField )
        % Written by LW
        % assumes 1 ms binned data in R struct
            pre_mo_ms = pre_post_mo_times( 1 );
            post_mo_ms = pre_post_mo_times( 2 );
            
            % trial length
            n_samples = numel( -1*pre_mo_ms:post_mo_ms ) - 1;
            
            % number of tot trials
            n_tot_trials = size( keep_trials, 2 );

            % number of keep trials
            n_keep_trials = sum( keep_trials );

            % find indices of keep trials
            keep_trial_inds = find( keep_trials == 1 );
            
            % number of channels
            n_channels = size( r.r(1).( fieldToGen ), 1 );

            % intialize logical vector to index spikes
            %keep_channels = true( 1, tot_channels );
            % set exclude channels for removal
            %keep_channels( exclude_channels ) = false;
            
            % initialize output counts matrix
            counts = zeros( n_keep_trials, n_channels, n_samples );
            % intialize output condition Id vector
            %condition_id = zeros( 1, n_keep_trials );
            % intialize output timeVecMs
            tmp = 1:( numel( -1*pre_mo_ms:post_mo_ms ) - 1 );
            time_vec_ms = tmp - pre_mo_ms;
            
            for i=1:n_keep_trials
                %i_trial = keep_trial_inds( i );
                i_trial = i;
                % get spikes
                data = r.r( i_trial ).( fieldToGen );
                %trial_cond = r.r( i_trial ).conditionId;
                % get alignment idx
                align_idx = r.r( i_trial ).( alignField );
                % get start and end idx for trials
                start_idx = align_idx - pre_mo_ms;
                end_idx = align_idx + post_mo_ms - 1;
                counts( i, :, : ) = data( :, start_idx:end_idx );
                %condition_id( i ) = trial_cond;
            end
            data = struct();
            data.counts = counts;
            %data.timeVecMs = time_vec_ms;
            %data.conditionId = condition_id;
        end

        function [ data ] = generate_lfads_data( r, pre_post_mo_times, keep_trials, exclude_channels )
        % Written by LW
        % assumes 1 ms binned data in R struct
            pre_mo_ms = pre_post_mo_times( 1 );
            post_mo_ms = pre_post_mo_times( 2 );
            
            % trial length
            n_samples = numel( -1*pre_mo_ms:post_mo_ms ) - 1;
            
            % number of tot trials
            n_tot_trials = size( keep_trials, 2 );

            % number of keep trials
            n_keep_trials = sum( keep_trials );

            % find indices of keep trials
            keep_trial_inds = find( keep_trials == 1 );
            
            % original number of channels
            tot_channels = size( r.r(1).spikes, 1 );
            % number of channels after high xcorr removal
            n_channels = tot_channels - numel( exclude_channels );

            % intialize logical vector to index spikes
            keep_channels = true( 1, tot_channels );
            % set exclude channels for removal
            keep_channels( exclude_channels ) = false;
            
            % initialize output counts matrix
            counts = zeros( n_keep_trials, n_channels, n_samples );
            % intialize output condition Id vector
            condition_id = zeros( 1, n_keep_trials );
            % intialize output timeVecMs
            tmp = 1:( numel( -1*pre_mo_ms:post_mo_ms ) - 1 );
            time_vec_ms = tmp - pre_mo_ms;
            
            for i=1:n_keep_trials
                i_trial = keep_trial_inds( i );
                % get spikes
                spikes = r.r( i_trial ).spikes;
                trial_cond = r.r( i_trial ).conditionId;
                % get movement onset idx
                move_on_idx = r.r( i_trial ).moveOnset;
                % get start and end idx for trials
                start_idx = move_on_idx - pre_mo_ms;
                end_idx = move_on_idx + post_mo_ms - 1;
                counts( i, :, : ) = spikes( keep_channels, start_idx:end_idx );
                condition_id( i ) = trial_cond;
            end
            data = struct();
            data.counts = counts;
            data.timeVecMs = time_vec_ms;
            data.conditionId = condition_id;
        end

        function test = findMovementOnsetTimes( r, speedThresh, goCueSpec )
        % DEPRECATED
        % please use find_move_onset_times
            disp( 'WARNING: This function is deprecated. Use find_move_onset_times for any future uses')
            nTrials = size( r.r, 2 );
            figure()
            for itrial = 1:nTrials
                vx = r.r( itrial ).kin( 3, : );
                vy = r.r( itrial ).kin( 4, : );
                %Calculate speed of monkey arm
                s = sqrt( vx.^2 + vy.^2 );
                % We only care about the changes in speed that occur after
                % the go cue, so we'll set all preceding indices to zero.
                switch goCueSpec
                  case 'rel'
                    goCue = r.r( itrial ).goCue;
                    s( 1, 1:goCue ) = 0;
                  case 'abs'
                    absGoCue = r.r( itrial ).goCue;
                    absTargOn = r.r( itrial ).targOn;
                    goCue = absGoCue - absTargOn;
                    test(itrial) = goCue; % DEBUG
                    %s( 1, 1:goCue ) = 0;
                  otherwise
                    disp('ERROR: Could not read goCueSpec option.')
                end
                % find within trial maximum speed
                maxSpeed = max( s );
                % find normalized speed for trial
                normSpeed = s / maxSpeed;
                plot( normSpeed )
                hold on
                % calculate index of when trial speed exceeds movement onset threshold
                moveOnsetIdx = find( normSpeed >= speedThresh, 1 );
                r.r( itrial ).moveOnset = moveOnsetIdx;
            end %for
        end %findMovementOnsetTimes
        function [ ] = plotTrialSpeed( r, trialNum, sampleRate, normFlag )
        % DEPRECATED
        % please use plot_speed
            disp( 'WARNING: This function is deprecated. Use plot_speed for any future uses')

            figure()
            nTrials = numel( trialNum );
            for itrial = 1:nTrials
                % get velocities 
                vX = r.r( trialNum( itrial) ).kin( 3, : );
                vY = r.r( trialNum( itrial) ).kin( 4, : );

                % calculate speed
                speed = sqrt( vX.^2 + vY.^2 );
                

                % find trial dataLength                
                trialDataLength = size( speed , 2 );
                
                % compare datalength of trials and return longest
                if itrial > 1
                    if trialDataLength > dataLength
                        dataLength = trialDataLength;
                    end
                else
                    dataLength = trialDataLength;
                end
                
                % if user specifies to normalize, normalize speed
                if normFlag
                    maxSpeed = max( speed );
                    speed = speed / maxSpeed;
                    plotTitle = 'Normalized Speed';
                    yLabel = 'Normalized Speed';
                else
                    plotTitle = 'Speed';
                    yLabel = 'Speed';
                end
                plot( speed )
                hold on
            end % for itrial
            
            % plot options
            disp( 'Plotting trials:' )
            if nTrials > 10
                disp( 'That is a lot of trials...')
            else
                disp( trialNum )
            end
            
            time = ( [ 1:dataLength ] * 1000 ) ./ sampleRate;
            title( plotTitle )
            xTickVec = 1:1000:dataLength;
            set(gca,'XTick',xTickVec);
            
            numTicks = numel( xTickVec );
            xTickVecLabels = cell( 1, numTicks );
            for i = 1:numTicks
                tickTime = time( xTickVec( i ) );
                xTickVecLabels{ i } = num2str( round( tickTime ) );
            end
            set( gca, 'XTickLabels', xTickVecLabels );
            xlabel( 'Time (ms)' )
            ylabel( yLabel )
            
        end % plotTrialSpeed
                function [ data ] = add_output_from_lfads( r, run, field_name_prefix )

            run.loadSequenceData();
            run.loadPosteriorMeans();
            run.addPosteriorMeansToSeq();
            
            % get sequence data
            seq_struct = run.sequenceData{ 1 };

            % number of trials
            n_trials = numel( seq_struct );
            
            % number of factors
            n_factors = size( seq_struct( 1 ).factors, 1 );
            
            % number of rates
            n_rates = size( seq_struct( 1 ).rates, 1 );
            
            % number of samples 
            n_samples = size( seq_struct( 1 ).rates, 2 );
            
            % initialize output for factors
            output_factors = zeros( n_trials, n_factors, n_samples );
            % intialize output for rates
            output_rates = zeros( n_trials, n_rates, n_samples );
            % intialize output for spiking data ( 1 ms bin )
            %output_spikes = zeros( n_trials, size( seq_struct( 1 ).y, 1 ), size( seq_struct( 1 ).y, 2 ) );

            % for each trial in seq struct
            for i=1:numel( seq_struct )
                % add factors and rates to 3-d matrix for output
                r.r( i ).( [ field_name_prefix, '_rates' ] ) = seq_struct( i ).rates;
                r.r( i ).( [ field_name_prefix, '_factors' ] ) = seq_struct( i ).factors;
            end
        end
    end %methods
end %classdef
%{
        function [ trial_conds, keep_trials ] = plot_move( r, pre_post_mo_ms , reach_time_max, max_speed_min, plot_title )
        % extract pre movement onset and post movement onset times
            pre_mo_ms = pre_post_mo_ms( 1 );
            post_mo_ms = pre_post_mo_ms( 2 );

            % determine number of trials
            n_trials = size( r.r, 2 );
            
            % intialize figure
            figure( 303 )
            
            % intialize count
            count = 1;
            
            % set x vector to be aligned to movement onset
            tmp = 1:( numel( -1*pre_mo_ms:post_mo_ms ) - 1 );
            plot_range = tmp - pre_mo_ms;

            % initialize color map
            cmap = hsv;
            numConds = numel( unique( [ r.r.conditionId ] ) );
            keep_trials = false( 1, numel( r.r ) ); 
            % for each trial
            for i = 1: n_trials
                % extract x and y velocities
                %v_x = r.r( i ).kin( 3, : );
                %v_y = r.r( i ).kin( 4, : );

                %s = sqrt( v_x.^2 + v_y.^2 );

                % calculate speed
                s = sqrt( sum( r.r(i).kin(3:4, :).^2 ) );
                
                % find value and index max speed
                [max_s, max_s_idx] = max( s );
                %max_s_idx = find( s == max_s );

                % get condition id
                icond = r.r( i ).conditionId;
                % get color for trace
                clrind = floor( (icond - 1) / numConds * size(cmap, 1) ) + 1;
                % get move onset
                move_onset_idx = r.r( i ).moveOnset;
                %                go_cue_idx = r.r( i ).relGoCue;
                % calculate reach time and move delay
                %                reach_time = max_s_idx - move_onset_idx;
                %                move_delay = move_onset_idx - go_cue_idx;
                % plot trial if trace meets conditions
                try
                    %                    if max_s >= max_speed_min
                        % plot speed trace
                        h = plot( plot_range, s( move_onset_idx - pre_mo_ms: move_onset_idx + post_mo_ms - 1) );
                        % add condition to keep trial condition vector
                        trial_conds( count ) = icond;
                        % change keep trial flag to true
                        keep_trials( i ) = true;
                        count = count + 1;
                        %                    end
                catch ME
                    disp( 'skipping trial' )
                end
                
                set( h, 'Color', cmap( clrind, : ) );
                %                    reach_time_keep( count ) = reach_time;
                %                    move_delay_keep( count ) = move_delay;
                hold on
            end
            ylabel( 'Speed (cm/s)')
            xlabel( 'Time Aligned to Move Onset (ms)' )
            title( plot_title )
            axis( 'tight' )
            %            figure( 202 )
            %            hist( reach_time_keep )
            %            axis( 'tight' )
            %            title( 'Time Delay between Move Onset and Max Speed' )
            %            figure( 101 )
            %            hist( move_delay_keep )
            %            axis( 'tight' )
            %            title( 'Time Delay between Go Cue and Move Onset' )
            disp( count )
        end
%}