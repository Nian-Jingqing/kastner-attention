function [index] = cycleIndex(origin_index)
%CYCLEINDEX say the indices could only be 1-4. If you have a index of 5,
%you returns a 1. If you have a index of 7, you returns a 3
%   Detailed explanation goes here
if rem(origin_index,4) == 0
    index = 4;
else
    index = rem(origin_index, 4);
end
end

