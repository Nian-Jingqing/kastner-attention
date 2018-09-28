classdef Continuous < handle
% class to hold / process continuous data
    properties
        data % struct, all fields should have the same first dimension (time)
        dtMS % sample rate of data, in MS
    end
    
    methods
        % constructor
        function c = Continuous( data, dtMS )
            if exist( 'data', 'var' )
                c.data = data;
            end
            if exist( 'dtMS', 'var' )
                c.dtMS = dtMS;
            end
        end
        
        function smoothField(c, infield, outfield, sigma)
        % apply a gaussian kernel to smooth the continuous data
            c.data.( outfield ) = Utils.gaussianSmooth( c.data.( infield ), sigma, 1, 0 );
        end
        
        function trials = makeTrialsFromData(c, start, stop, trialFields)
        % turn the continuous data into a struct where each element is a trial
        % start: indices into the continuous data where trials start
        % stop: indices where trials end
        % trialFields: any data to be stored with each trial
            assert( numel(start) == numel(stop), 'huh??');
            assert( all( ( stop - start ) > 0 ), 'huh??');
            if exist('trialFields', 'var')
                assert( numel(stop) == numel(trialFields), 'huh??');
            end
            
            % split all the fields from continuous data into a struct array
            cfields = fields( c.data );
            % also add all the trial fields
            tfields = fields(trialFields);

            trials = [];
            % iterate over all trials
            for itrial = 1 : numel( start )
                % first add the continuous data
                for ifield = 1 : numel ( cfields )
                    dat = c.data.( cfields{ ifield } )( ...
                        start(itrial) : stop(itrial)-1, : );
                    % for trial-ized data, we store with time as 2nd dimension
                    trials( itrial ).( cfields{ifield} ) = dat';
                end
                % next add the trial info
                for ifield = 1 : numel( tfields )
                    trials( itrial ).( tfields{ifield} ) = ...
                        trialFields( itrial ).( tfields{ifield} );
                end
                
            end
        end


    end


end
