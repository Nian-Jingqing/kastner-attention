function alignedStruct = alignedStructToSingleTrial_lfads( singleTrialStruct, fieldsToConvert )
% function Data = singleTrialStructToAlignedStruct( singleTrialStruct, fieldsToConvert )
%   singleTrialStruct should be a struct array
%    for every field we want to convert, its leading dimension should be an individual trial
%    those trials will be converted to individual conditiosn

if ~iscell( fieldsToConvert )
    fieldsToConvert = { fieldsToConvert };
end

% get all the fields in the struct
allFields = fields( singleTrialStruct );
% get the ones we arent "converting", those will just be stored in the new struct as-is
otherFields = setdiff( allFields, fieldsToConvert );

%% take an aligned struct, and convert it to the weird generic format that jPCA accepts
numConds = numel( singleTrialStruct );

outIndex = 0;
alignedStruct = struct;
% iterate over all conditions
for c = 1:numConds
    % get the number of trials for this condition. later, check that is constant for all fields we're converting
    numTrialsThisCondition = size( singleTrialStruct( c ).( fieldsToConvert{ 1 } ), 1 );

    % iterate over trials
    for ntrial = 1:numTrialsThisCondition
        outIndex = outIndex + 1;
        
        %iterate over fields being converted
        for nf = 1:numel( fieldsToConvert )
            origsize = size( singleTrialStruct( c ).( fieldsToConvert{ nf } ) );
            % check this field has the same number of trials as the first one
            assert( origsize(1) == numTrialsThisCondition, 'error: number of trials is not constant' );

            % store down this trial's data
            %  the reshape and squeeze are important because we don't know if the original data is e.g. 2-dimensional, 3-dimensional, etc
            alignedStruct( outIndex ).( fieldsToConvert{ nf } ) = ...
                reshape( squeeze( singleTrialStruct( c ).( fieldsToConvert{ nf } )( ntrial, : ) ), origsize( 2:end ) );
        end
        
        % store the remaining fields
        for nf = 1:numel( otherFields )
            alignedStruct( outIndex ).( otherFields{ nf } ) = singleTrialStruct( c ).( otherFields{ nf } );
        end % otherFields
    end  %numTrialsThisCondition
    
end % numConds
