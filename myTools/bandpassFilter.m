function  R = bandpassFilter( R, fieldToFilter, outputField, filtHighCutoff, filtLowCutoff, Fs )
%BANDPASSFILTER Summary of this function goes here
%   Detailed explanation goes here
    
    % need to filter the LFP, but we need it to be continuous.
    LFP = [R.(fieldToFilter)];
    % filter the LFP
    [b,a] = butter(4, [filtLowCutoff filtHighCutoff] / (Fs / 2), 'bandpass');
    LFP2 = filtfilt(b, a, LFP');

    % we didn't have defined LFP at the edges, or at any clock discontinuities
    % add 100 ms of nan buffer at those times
    LFP2(1:100,:) = nan;
    LFP2(end-(100:0),:) = nan;

    % now place the filtered LFP back into the struct
    nTimes = size(R(1).(fieldToFilter), 2);
    chopBlock = 1:nTimes;
    for trial = 1: length(R)
        R(trial).(outputField) = LFP2(chopBlock, :)';
        chopBlock = chopBlock + nTimes;
    end

%         discont = find(diff(clock)~=1);
% 
%         for nd = 1:numel(discont)
%             LFP2(discont(nd)+(-100:100), :) = nan;
%         end
% 
%         % now place the filtered LFP back into the struct
%         currind = 0;
%         for ntr = 1:numel(R)
%             tinds = 1:numel(R(ntr).clock);
%             R(ntr).LFP = LFP2(tinds + currind, :)';
%             currind = currind+tinds(end);
%         end
    
end

