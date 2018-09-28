function [matrixFiltered] = bandpassFilter_continuous(matrixToFilter, filtHighCutoff, filtLowCutoff, Fs)
%BANDPASSFILTER_CONTINUOUS Summary of this function goes here
%   Detailed explanation goes here
    % filtering
    [b,a] = butter(4, [filtLowCutoff filtHighCutoff] / (Fs / 2), 'bandpass');
    matrixFiltered = filtfilt(b, a, matrixToFilter');

    % we didn't have defined LFP at the edges, or at any clock discontinuities
    % add 100 ms of nan buffer at those times
    matrixFiltered(1:100,:) = nan;
    matrixFiltered(end-(100:0),:) = nan;
    matrixFiltered = matrixFiltered';

end

