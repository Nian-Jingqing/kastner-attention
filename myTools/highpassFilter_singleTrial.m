function  outputMatrix = highpassFilter_singleTrial( inputMatrix,  filtLowCutoff, Fs )
%BANDPASSFILTER Summary of this function goes here
%   Detailed explanation goes here
[b,a] = butter(4, filtLowCutoff / (Fs / 2), 'high');
LFP = inputMatrix;
addFront = repmat(LFP(:, 1), [1, 20]);
addBack = repmat(LFP(:,end), [1,20]);
LFP = [addFront, LFP, addBack];
LFP2 = filtfilt(b, a, LFP');
outputMatrix = LFP2(21: end-20, :)';