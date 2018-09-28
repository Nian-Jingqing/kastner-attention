classdef outputData < Movement.movementData
    methods
        % constructor
        function obj = outputData( varargin )
            obj = obj@Movement.movementData ( varargin{:} );
        end
    end %methods
end %classdef
