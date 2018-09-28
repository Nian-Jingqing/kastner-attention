classdef TrialType < uint32
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    enumeration
        InterTrial (0)
        Success (1)
        Failure (2)
    end
    methods
        function c = eq(a,b)
            c = isequal(uint32(a),uint32(b));
        end
    end
    
end

