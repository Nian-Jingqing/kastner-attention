classdef Rstruct < handle & matlab.mixin.Copyable
    properties
        r % struct, each element is a different trial
    end

    methods
        % constructor
        function r = Rstruct( data )
            if exist( 'data', 'var' )
                r.r = data;
            end
        end
        function [alignedStruct] = getAligned( r, conditionField, fieldsToAlign, prePostTimes, alignField )
        % extract all trials of a given condition, aligned by a certain field
        %    conditionField - name of a field in the struct to use to split into conditions
        %    fieldsToAlign - cell array of fields that you want to average
        %    prePostTimes - 2 element vector with time before and time after the align point
        %    alignField - if the data needs to be re-aligned, what field to use

            % quick checks
            assert( ischar( conditionField ), ...
                    'second argument is field to use to split trials into conditions' );
            % fieldsToAlign is a cell array
            if ~iscell( fieldsToAlign )
                fieldsToAlign = { fieldsToAlign };
            end

            % default for 'alignField' is empty
            if ~exist( 'alignField', 'var' )
                alignField = [];
            end

            % split trials by condition
            codes = [ r.r.( conditionField ) ];
            [ucodes, ~, trialInds] = unique(codes);

            % how long is the sequence
            sequenceLength = prePostTimes(1) + prePostTimes(2);
            sequenceOffsets = ( -prePostTimes(1) : prePostTimes(2)-1 ) + 1;

            % get the data dimensionality for each field
            for ifield = 1:numel( fieldsToAlign )
                dataDim(ifield) = size( r.r(1).( fieldsToAlign{ifield} ), 1 );
            end % ifield

            % iterate over conditions
            for icond = 1:numel( ucodes )
                % get indices of trials for this condition
                thesetrials = find( trialInds == icond );

                % store some info about the condition
                alignedStruct(icond).conditionCode = ucodes(icond);
                alignedStruct(icond).numTrials = numel(thesetrials);
                alignedStruct(icond).origTrialIdx = thesetrials;
                
                % get all the trials for each field
                for ifield = 1:numel( fieldsToAlign )
                    fname = fieldsToAlign{ ifield };

                    % place to store all trials
                    alignedStruct( icond ).( fname ) = zeros( numel(thesetrials), dataDim(ifield), ...
                                                              sequenceLength );
                    % iterate over trials for this condition
                    for itrial = 1:numel(thesetrials)
                        trialnum = thesetrials( itrial );
                        % get the alignpoint for this trial
                        if isempty( alignField )
                            alignPoint = 1;
                        else
                            alignPoint = r.r( trialnum ).( alignField );
                        end

                        % what indices do we get for this trial?
                        indsToGet = alignPoint + sequenceOffsets;
                        assert( all( indsToGet > 0 ), 'error in indexing - some negative values' );
                        % extract the data
                        data = r.r( trialnum ).(fname)( : , indsToGet );
                        % put the data into the allTrials array
                        alignedStruct( icond ).( fname )(itrial, :, : ) = data;

                    end % itrial

                end % ifield

                % store the time vector
                alignedStruct( icond ).timeVecMs = sequenceOffsets;
            end % icond

        end % getAligned

        function [avgStruct, stdStruct] = alignAndAverage( r, conditionField, fieldsToAverage, prePostTimes, alignField )
        %  average across all trials of a given condition (e.g., to create a PSTH)
        %    conditionField - name of a field in the struct to use to split into conditions
        %    fieldsToAverage - cell array of fields that you want to average
        %    prePostTimes - 2 element vector with time before and time after the align point
        %    alignField - if the data needs to be re-aligned, what field to use
            % quick checks
            assert( ischar( conditionField ), ...
                    'second argument is field to use to split trials into conditions' );
            % fieldsToAverage is a cell array
            if ~iscell( fieldsToAverage )
                fieldsToAverage = { fieldsToAverage };
            end

            % default for 'alignField' is empty
            if ~exist( 'alignField', 'var' )
                alignField = [];
            end

            % how long is the sequence
            sequenceLength = prePostTimes(1) + prePostTimes(2);
            sequenceOffsets = ( -prePostTimes(1) : prePostTimes(2)-1 ) + 1;

            % call anther function to get the aligned data
            alignedStruct = r.getAligned( conditionField, fieldsToAverage, prePostTimes, alignField );

            % do the averaging below...

            % iterate over conditions
            for icond = 1:numel( alignedStruct )
                % store the condition info down
                avgStruct( icond ).conditionCode = alignedStruct(icond).conditionCode;
                avgStruct( icond ).numTrials = alignedStruct(icond).numTrials;
                stdStruct( icond ).conditionCode = alignedStruct(icond).conditionCode;
                stdStruct( icond ).numTrials = alignedStruct(icond).numTrials;

                % average / std the values for each field
                for ifield = 1:numel( fieldsToAverage )
                    fname = fieldsToAverage{ ifield };

                    avgStruct( icond ).( fieldsToAverage{ifield} ) = ...
                        squeeze( mean( alignedStruct( icond ).( fieldsToAverage{ifield} )  ) );
                    stdStruct( icond ).( fieldsToAverage{ifield} ) = ...
                        squeeze( std( alignedStruct( icond ).( fieldsToAverage{ifield} )  ) );
                end % ifield

                % store time vector, if it exists
                if isfield( alignedStruct, 'timeVecMs' )
                    avgStruct( icond ).timeVecMs = alignedStruct.timeVecMs;
                    stdStruct( icond ).timeVecMs = alignedStruct.timeVecMs;
                end
            end % icond

        end % alignAndAverage

        function [avgStruct, stdStruct] = alignAndAverage_lfads( r, conditionField, fieldsToAverage, prePostTimes, lfadsTimes, alignField )
        %  average across all trials of a given condition specifically with lfads output (e.g., to create a PSTH)
        %    conditionField - name of a field in the struct to use to split into conditions
        %    fieldsToAverage - cell array of fields that you want to average
        %    prePostTimes - 2 element vector with time before and time after the align point
        %    lfadsTimes- a vector of the times corresponding to the lfads output bins
        %    alignField - if the data needs to be re-aligned, what field to use
            % quick checks
            assert( ischar( conditionField ), ...
                    'second argument is field to use to split trials into conditions' );
            % fieldsToAverage is a cell array
            if ~iscell( fieldsToAverage )
                fieldsToAverage = { fieldsToAverage };
            end

            % default for 'alignField' is empty
            if ~exist( 'alignField', 'var' )
                alignField = [];
            end

            % how long is the sequence
            sequenceLength = prePostTimes(1) + prePostTimes(2);
            sequenceOffsets = ( -prePostTimes(1) : prePostTimes(2)-1 ) + 1;

            % call anther function to get the aligned data
            alignedStruct = r.getAligned( conditionField, fieldsToAverage, prePostTimes, alignField );

            % do the averaging below...

            % iterate over conditions
            for icond = 1:numel( alignedStruct )
                % store the condition info down
                avgStruct( icond ).conditionCode = alignedStruct(icond).conditionCode;
                avgStruct( icond ).numTrials = alignedStruct(icond).numTrials;
                stdStruct( icond ).conditionCode = alignedStruct(icond).conditionCode;
                stdStruct( icond ).numTrials = alignedStruct(icond).numTrials;

                % average / std the values for each field
                for ifield = 1:numel( fieldsToAverage )
                    fname = fieldsToAverage{ ifield };

                    avgStruct( icond ).( fieldsToAverage{ifield} ) = ...
                        squeeze( mean( alignedStruct( icond ).( fieldsToAverage{ifield} )  ) );
                    stdStruct( icond ).( fieldsToAverage{ifield} ) = ...
                        squeeze( std( alignedStruct( icond ).( fieldsToAverage{ifield} )  ) );
                end % ifield

                % store time vector, if it exists
                if isfield( alignedStruct, 'timeVecMs' )
                    %We want to overwrite the times to be the lfadsTimes
                    avgStruct( icond ).timeVecMs = lfadsTimes;
                    stdStruct( icond ).timeVecMs = lfadsTimes;
                end
            end % icond

        end % alignAndAverage_lfads

        function rbinned = binData( r, fieldsToBin, binSize )
        % iterate over all trials
            for itrial = 1:numel( r.r )

                % iterate over all the fields we want to bin
                for nf = 1:numel( fieldsToBin )
                    d = r.r( itrial ).( fieldsToBin{ nf } );
                    dataDim = size( d, 1);

                    % trim any data that won't fit in an even number of bins
                    dataLength = size( d, 2);
                    numBins = floor( dataLength / binSize(nf) );
                    pointsToKeep = numBins * binSize(nf);
                    d = d( :, 1:pointsToKeep );

                    % reshape, sum, and squeeze
                    d2 = reshape( full( d ), dataDim, binSize(nf), numBins );
                    dataBinned = squeeze( sum( d2, 2 ) );

                    % if the original data is 1-dimensional, there are some problems with squeeze
                    %  correction for that is below
                    if dataDim == 1
                        dataBinned = dataBinned(:)';
                    end

                    % assign that data
                    rbinned ( itrial ).( fieldsToBin{ nf } ) = dataBinned;

                end % nf
            end % itrial
        end % binData

        function [alignedStruct] = getAlignedTrials( r, fieldsToAlign, prePostTimes, alignField )
        % extract all trials of a given condition, aligned by a certain field
        %    fieldsToAlign - cell array of fields that you want to average
        %    prePostTimes - 2 element vector with time before and time after the align point
        %    alignField - if the data needs to be re-aligned, what field to use

            % fieldsToAlign is a cell array
            if ~iscell( fieldsToAlign )
                fieldsToAlign = { fieldsToAlign };
            end

            % default for 'alignField' is empty
            if ~exist( 'alignField', 'var' )
                alignField = [];
            end

            % how long is the sequence
            sequenceLength = prePostTimes(1) + prePostTimes(2);
            sequenceOffsets = ( -prePostTimes(1) : prePostTimes(2)-1 ) + 1;

            % get the data dimensionality for each field
            for ifield = 1:numel( fieldsToAlign )
                dataDim(ifield) = size( r.r(1).( fieldsToAlign{ifield} ), 1 ); % dimensionality of the data
            end % ifield

            % get the total number of trials
            numTrials = numel(r.r);

            % iterate over trials
            for itrial = 1:numTrials
                % get alignpoint for this trial
                if isempty( alignField )
                    alignPoint = 1;
                else
                    alignPoint = r.r( itrial ).( alignField );
                end
                indsToGet = alignPoint + sequenceOffsets;
                assert( all( indsToGet > 0), 'error in indexing - some negative values');

                % iterate over fields for each trial
                for ifield = 1:numel( fieldsToAlign )
                    fname = fieldsToAlign{ ifield };
                    % place to store all data for a particular field in a trial
                    alignedStruct(itrial).( fname ) = zeros( dataDim(ifield), sequenceLength );
                    % extract the data
                    data = r.r( itrial ).(fname)( :, indsToGet );
                    % put the data into the allTrials array
                    alignedStruct( itrial ).( fname ) = data;

                    r.r(itrial).(fname) = data;


                end % ifield
            end % itrials
        end % getAlignedTrials

        
        function smoothFieldInR( r, infield, outfield, sigma, rebinSize )
        % apply a gaussian kernel to smooth the continuous data
           for itrial = 1:numel(r.r)
               r.r( itrial ).(outfield) = (Utils.gaussianSmooth(r.r( itrial ).(infield)', sigma, rebinSize,0))';
           end % itrial
        end % smoothFieldInR
        
        function postSmoothCutOff( r, infield, outfield, cutOff)
            for itrial = 1:numel(r.r)
                r.r( itrial ).(outfield) = r.r( itrial ).( infield )(:,(cutOff+1):(end-cutOff));
            end % itrial
        end % postSmoothCutOff

    end % methods
end % classdef
