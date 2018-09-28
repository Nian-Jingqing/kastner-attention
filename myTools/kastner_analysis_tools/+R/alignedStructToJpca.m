function Data = alignedStructToJpca( alignedStruct, varargin )
% function Data = alignedStructToJpca( alignedStruct, neuralField, conditionField )
%   alignedStruct should be a struct array, one element per "condition"

% use matlab's "inputParser" to handle input to this function
p = inputParser();
% neuralField must be passed in as a 2nd argument
p.addRequired( 'neuralField', @ischar );
% conditionField is an optional named argument. [i.e., (..., 'conditionField', 'thisFieldHasMyConditionInfo' );
%   default is []
p.addParameter( 'conditionField', [], @ischar );
% timesField is an optional named argument. [i.e., (..., 'timesField', 'thisFieldHasMyTimeInfo' );
%   default is 'timeVecMs'
p.addParameter( 'timesField', 'timeVecMs', @ischar );

p.parse( varargin{:} );
res = p.Results;

%% take an aligned struct, and convert it to the weird generic format that jPCA accepts
numConds = numel( alignedStruct );
numNeurons = size( alignedStruct(1).( res.neuralField ), 1 );
numTimes = size( alignedStruct(1).( res.neuralField ), 2 );
% size of alignedStruct.(neuralField) should be the same as the time vector
assert( numTimes == numel( alignedStruct(1).( res.timesField ) ), 'error: conflicting ideas re: the time dimension')


for c = 1:numConds
    % time vector
    Data( c ).times = alignedStruct( c ).( res.timesField );
    % neural data
    Data( c ).A = alignedStruct( c).( res.neuralField )';

    % store condition code if it was passed in
    if ~isempty( res.conditionField )
        Data( c ).( res.conditionField ) = alignedStruct( c ).( res.conditionField );
    end
end
