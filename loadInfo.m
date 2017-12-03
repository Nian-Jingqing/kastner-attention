function [rt, cueLoc, window] = loadInfo(filename)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
S = load(filename);
window = S.targetDim.window;
rt = S.UE.rt(S.UE.isHoldTrial);%need to figure out the rt for the right trials
cueLoc = S.UE.cueLoc(S.UE.isHoldTrial);
%need to figure out the cue location for the right trials
end

