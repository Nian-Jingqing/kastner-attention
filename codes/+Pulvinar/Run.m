classdef Run < LFADS.Run
    methods
        function r = Run(varargin)
           r@LFADS.Run(varargin{:});
        end

        function tf = usesDifferentDataForAlignment(r)  %#ok<MANU>
            % tf = usesDifferentDataForAlignment()
            %
            % Returns true if you would like the Run to call your
            % convertDatasetToSequenceStruct with a mode == 'alignment'
            % argument when constructing the alignment matrices. This
            % allows you to specify a different set of trials used for
            % constructing alignment matrices, e.g. only correct trials.

                tf = true; %uncomment this if using different datasets for alignment
                   %      tf = false;
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
            
            % get run params
            
            
            
            
         %%%%%%%%%%%% for using a different dataset for alignment matrix %%%%%%%%%
                          par = r.params;
                    if strcmp(mode, 'export')
                        data = dataset.loadData();
                        data = data.combinedData;
                
                       trial_time_ms = 500;
                       trial_olap_ms = 100;
                       out = data.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );
                 
                   elseif strcmp(mode, 'alignment')
                       data = dataset.loadData();
                       data = data.combinedData;
                       spikeTensor = cat(3, data.R.spikeCounts);
                       if ~isempty( par.nIndices )
                           spikeTensor = spikeTensor( par.nIndices, : , : );
                       end
                       out.counts = permute( spikeTensor, [3 1 2] ); % use "permute" to make the trials to 
                       %the first dimension. Now the out.counts is a nTrials x
                       %nChannels x nTime matrix
                        %       out.timeVecMs = 1: size(out.counts,3);
                        %       out.conditionId = (ones(1, size(out.counts,1)))';
                       out.conditionId = ([data.R.condition])';
                   end % must have this end
         
                end %must have this end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%############ runs for pre-aligned data but with external inputs ##################

        %par = r.params;
   
        %data = dataset.loadData();

        %for itrial = 1 : numel(data.R)
        %    data.R(itrial).spikeCounts = full(data.R(itrial).spikeCounts);
        %    data.R(itrial).externalInputs = full(data.R(itrial).externalInputs);
        %end
        
        %spikeTensor = cat(3, data.R.spikeCounts); % concatenate R struct to 3D matrix
        %EITensor = cat(3, data.R.externalInputs);
        %% if user specifies neurons, extract those indices 
        %if ~isempty( par.nIndices )
        %    spikeTensor = spikeTensor( par.nIndices, : , : );
        %    EITensor = EITensor( par.nIndices, :, : );
        %end
        %% Note, the matrix is nChannels x nTime x nTrials
        %out.counts = permute( spikeTensor, [3 1 2] ); % use "permute" to make the trials to
        %out.externalInputs = permute( EITensor, [3 1 2] );
        %%the first dimension. Now the out.counts is a nTrials x
        %%nChannels x nTime matrix
        %   %     out.timeVecMs = 1: size(out.counts,3);
        %    %      out.conditionId = (ones(1, size(out.counts,1)))';
        %out.conditionId = ([data.R.condition])';
        %   %      out.externalInputs = [];
        %% out.counts = data.spikes; % for lorenz dataset          
        %% out.timeVecMs = data.timeMs; % for lorenz dataset
        %% out.conditionId = data.conditionId; % for lorenz dataset
        %% out.truth = data.true_rates; % for lorenz dataset
        %end
            
            
            
           %%%%%%%%%%%%%%%%%%%%%%%% runs for no external inputs %%%%%%%%%%%%%%%%%
        
           %par = r.params;
           %
           %data = dataset.loadData();
        
           %spikeTensor = cat(3, data.R.spikeCounts); % concatenate R struct to 3D matrix. 
           %% if user specifies neurons, extract those indices 
           %if ~isempty( par.nIndices )
           %    spikeTensor = spikeTensor( par.nIndices, : , : );
           %end
           %% Note, the matrix is nChannels x nTime x nTrials
           %out.counts = permute( spikeTensor, [3 1 2] ); % use "permute" to make the trials to 
           %%the first dimension. Now the out.counts is a nTrials x
           %%nChannels x nTime matrix
           %     %     out.timeVecMs = 1: size(out.counts,3);
           %     %      out.conditionId = (ones(1, size(out.counts,1)))';
           % out.conditionId = ([data.R.condition])';
           %    %      out.externalInputs = [];
           %% out.counts = data.spikes; % for lorenz dataset          
           %% out.timeVecMs = data.timeMs; % for lorenz dataset
           %% out.conditionId = data.conditionId; % for lorenz dataset
           %% out.truth = data.true_rates; % for lorenz dataset
           %end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            
            
            
            
    end
end
