function  R = bandpassFilter_trialized( R, fieldToFilter, outputField, filtHighCutoff, filtLowCutoff, Fs )
%BANDPASSFILTER Summary of this function goes here
%   Detailed explanation goes here
    [b,a] = butter(4, [filtLowCutoff filtHighCutoff] / (Fs / 2), 'bandpass');
    for trial = 1:length(R)
        LFP = R(trial).(fieldToFilter);
        addFront = repmat(LFP(:, 1), [1, 20]);
        addBack = repmat(LFP(:,end), [1,20]);
        LFP = [addFront, LFP, addBack];
        LFP2 = filtfilt(b, a, LFP');
        R(trial).(outputField) = LFP2(21: end-20, :)';
    end
    
    %    % need to filter the LFP, but we need it to be continuous.
    % LFP = [R.(fieldToFilter)];
    % filter the LFP
    % [b,a] = butter(4, [filtLowCutoff filtHighCutoff] / (Fs / 2), 'bandpass');
    % LFP2 = filtfilt(b, a, LFP');

    % % we didn't have defined LFP at the edges, or at any clock discontinuities
    % % add 100 ms of nan buffer at those times
    % LFP2(1:100,:) = nan;
    % LFP2(end-(100:0),:) = nan;

    % % now place the filtered LFP back into the struct
    % trialStart = 1;
    % for trial = 1: length(R)
    %     trialEnd = trialStart + size(R(trial).(fieldToFilter),2) -1;
    %     R(trial).(outputField) = LFP2(trialStart:trialEnd, :)';
    %     trialStart = trialEnd + 1;
    % end
    
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


