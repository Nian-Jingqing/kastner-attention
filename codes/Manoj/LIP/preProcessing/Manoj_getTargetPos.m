function [targetPos] = Manoj_getTargetPos(cueType)
% This function compute targetPosition, based on the cued location
if ismember(cueType, [4, 8])
    targetPos = 4;
else
    targetPos = mod(cueType,4);
end
end