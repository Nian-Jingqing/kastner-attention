function [ rfLoc, spikeTimes ] = loadNeuronForContinuous( filename )
% For continuous
%   Detailed explanation goes here
S = load(filename);
spikeTimes = S.spikeTs*1000;
rfLoc = zeros(1, 2);
rfLoc(1) = S.inRFLoc;
rfLoc(2) = S.exRFLoc;
end

