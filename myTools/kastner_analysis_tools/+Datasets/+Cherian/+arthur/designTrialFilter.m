function [filt_main] = designTrialFilter(nbins,trialTime,trialOlap,filterType)
    % Design ascending linear filter
    leftLinFilt = 1:trialOlap;
    leftLinFilt = leftLinFilt/trialOlap;
    % Design descending linear filter
    rightLinFilt = trialOlap-1:-1:0;
    rightLinFilt = rightLinFilt/trialOlap;
        
    switch filterType
      case 0
        % Filter for first trial
        filt_main = ones(1,trialTime);
        filt_main(end-trialOlap+1:end) = rightLinFilt;
      case 1
        % Filter for middle trials                
        filt_main = ones(1,trialTime);
        filt_main(1:trialOlap) = leftLinFilt;
        filt_main(end-trialOlap+1:end) = rightLinFilt;
      case 2
        % Filter for last trial
        trialNewData = mod(nbins-trialTime,trialTime-trialOlap);
        nTrialIgnore = trialTime - trialNewData - trialOlap;
        trialIgnore = 1:nTrialIgnore;
        trialFilter = nTrialIgnore:nTrialIgnore+trialOlap-1;
        filt_main = ones(1,trialTime);
        filt_main(trialIgnore) = 0;
        filt_main(trialFilter) = leftLinFilt;
    end
end
