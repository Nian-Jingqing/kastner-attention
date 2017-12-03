function [spikeTimes, loc] = loadNeuron(filename)
%This function extracts the spike trains for all trials, receptive field, and the opposite
%RF for the neuron
%   Detailed explanation goes here
S = load(filename);
spikeTimes = S.targetDim.spikeTimes;
loc(1) = S.inRFLoc;
loc(2) = S.exRFLoc;
end

