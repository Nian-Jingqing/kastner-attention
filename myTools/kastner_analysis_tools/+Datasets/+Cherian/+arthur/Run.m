classdef Run < LFADS.Run
    methods
        function r = Run(varargin)
           r@LFADS.Run(varargin{:});
        end

        function out = generateCountsForDataset(r, dataset, mode, varargin) %#ok<INUSL,INUSD>
            % Generate binned spike count tensor for a single dataset.
            %
            % Parameters
            % ------------
            % dataset : :ref:`LFADS_Dataset`
            %   The :ref:`LFADS_Dataset` instance from which data were loaded
            %
            % mode (string) : typically 'export' indicating sequence struct
            %   will be exported for LFADS, or 'alignment' indicating that this
            %   struct will be used to generate alignment matrices. You can
            %   include a different subset of the data (or different time
            %   windows) for the alignment process separately from the actual
            %   data exported to LFADS, or return the same for both. Alignment
            %   is only relevant for multi-dataset models. If you wish to use
            %   separate data for alignment, override the method usesDifferentDataForAlignment
            %   to return true as well.
            %
            % Returns
            % ----------
            % out: a scalar struct with the following fields:
            %
            % .counts : nTrials x nChannels x nTime tensor
            %   spike counts in time bins in trials x channels x time. These
            %   should be total counts, not normalized rates, as they will be
            %   added during rebinning.
            %
            % .timeVecMs: nTime x 1 vector
            %   of timepoints in milliseconds associated with each time bin. You can start this
            %   wherever you like, but timeVecMs(2) - timeVecMs(1) will be
            %   treated as the spike bin width used when the data are later
            %   rebinned to match run.params.spikeBinMs
            %
            % .conditionId: nTrials x 1 vector
            %   of unique conditionIds. Can be cell array of strings or
        %   vector of unique integers.
        % get params
            par = r.params;
            % get i_blocks for run
            block = par.i_block;
            

            % get neurons to exclude from run
            exclude_neuron_path = dataset.exclude_neuron_path;
            
            % load file containing indices of neurons to remove 'all_neurons_to_remove'
            load( exclude_neuron_path )

            % get dataset type: 'center-out' or 'random-walk' ( for dataset loading purposes )
            dataset_type = par.dataset_type;
            
            % method of how trials are generated: 'align-chop' (typically 'center-out') or 'overlap-chop' (typically 'random-walk')
            trialize_method = par.trialize_method;

            run_type = par.run_type;

            % for each block in run
            for i = 1:numel( block )
                % get block number
                i_block = block( i );
                % load dataset
                dat = dataset.loadData( i_block, dataset_type );
                


                if strcmp( trialize_method, 'align-chop' ) % if trials generated using 'align-chop'
                    ap_threshold = par.ap_threshold;
                    pre_ap_ms = par.pre_ap_ms;
                    post_ap_ms = par.post_ap_ms;
                    ap_max_speed_min = par.ap_max_speed_min;
                    pre_post_ap_ms = [ pre_ap_ms post_ap_ms ];
                    % find alignment point for each trial
                    dat.find_align_point_times( ap_threshold, 'move_align' );
                    % determine which trials to keep for run
                    keep_trials = dat.find_align_chop_keep_trials( ap_threshold, ap_max_speed_min, pre_post_ap_ms, 'move_align' );
                    % generate LFADS counts data struct
                    data = dat.generate_align_chop_lfads_data( pre_post_ap_ms, keep_trials, run_type, all_neurons_to_remove, 'move_align' );
                    
                elseif strcmp( trialize_method, 'overlap-chop' ) % if trials generated using 'overlap-chop'
                    trial_time_ms = par.trial_time_ms;
                    trial_olap_ms = par.trial_olap_ms;
                    % determine which trials to keep for run (currently all)
                    keep_trials = dat.find_overlap_chop_keep_trials();
                    % generate LFADS counts data struct
                    data = dat.generate_overlap_chop_lfads_data( keep_trials, trial_time_ms, trial_olap_ms, run_type, all_neurons_to_remove );
                    
                else % if neither of the above options
                    assert( strcmp( trialize_method, 'align-chop' ) || strcmp( trialize_method, 'overlap-chop' ), ...
                            'ERROR: trialize_method must either be "align-chop" or "overlap-chop"' )
                end
                if i == 1
                    out.counts = data.counts;
                    out.timeVecMs = data.timeVecMs;
                    out.conditionId = data.conditionId;
                %out.truth = data.true_rates;
                else
                    out.counts = cat( 1, out.counts, data.counts );
                    out.conditionId = cat( 2, out.conditionId, data.conditionId );
                end
            end

        end
        
        function [ binned_data ] = rebin_data( r, out, binSize )
        % iterate over all trials
            for itrial = 1:size( out.counts, 1 );

                % iterate over all the fields we want to bin
                d = squeeze( out.counts( itrial, :, : ) );
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

                binned_data( itrial, :, : ) = dataBinned;
            end % itrial
        end % binData
        function [ factors, rates, condition_id ] = get_output_from_lfads( r, out_name )

            r.loadSequenceData();
            r.loadPosteriorMeans();
            r.addPosteriorMeansToSeq();
            
            % get sequence data
            seq_struct = r.sequenceData{ 1 };

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
                output_rates( i, :, : ) = seq_struct( i ).rates;
                output_factors( i, :, : ) = seq_struct( i ).factors;
                %output_spikes( i, :, : ) = seq_struct( i ).y;
            end

            % get vector of condition ids
            condition_id = [ seq_struct.conditionId ];

            out_folder = r.pathLFADSOutput;

            out_file = [ out_name '_output.mat' ];
            out_path = fullfile( out_folder, out_file );
            
            % rename factors and rates for output
            %factors = r.reshape_data_to_2d( output_factors, size( output_factors ) );
            %rates = r.reshape_data_to_2d( output_rates, size( output_rates ) );
            factors = output_factors;
            rates = output_rates;

            % rename 1ms binned spikes
            %binned_spikes_1ms = output_spikes;

            sample_rate_data = 1000 / r.params.spikeBinMs;
            disp( sprintf( 'writing file to %s\n', out_path ) );
            save( out_path, 'factors', 'rates', 'condition_id', 'sample_rate_data' );
        end
        function [ reshape_data ] = reshape_data_to_2d( r, data, data_shape )
        % get number of trials, dimensions, and samples for reshape
            n_trials = data_shape( 1 );
            n_dim = data_shape( 2 );
            n_samples = data_shape( 3 );
            
            % permute data for reshape
            perm_data = permute( data, [ 2 3 1 ] );
            
            % rehshape data and return
            reshape_data = reshape( perm_data, [ n_dim, n_trials*n_samples ] );
        end
    end
end
