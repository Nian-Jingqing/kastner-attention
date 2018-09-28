function [word_idx] = extractWordIndices(words,target_times,numSuccessTargets)
% Function to return indices to access data associated with:
%   (1) End of previous trial
%   (2) Start of trial
%   (3) Target 1 on
%   (4) Target 2 on
%   (5) Target 3 on
%   (6) End of trial
% in the words table to eventually extract relevant sample indices

    % Extract times from input word_idx vector
    t1 = target_times(1);
    t2 = target_times(2);
    t3 = target_times(3);
    % Extract word index for first target
    word_targ1_idx = find( words( :, 1 ) == t1 );
    % Determine trial start and end based off of target events
    word_start_idx = word_targ1_idx-1;
    assert(words(word_start_idx,2) == double(18),'ERROR: Start word not before first target.')
    word_end_idx = find(words(word_targ1_idx:end,2) == double(18),1,'first') + word_targ1_idx - 2;

    if numSuccessTargets >= 1
        word_targ2_idx = find(words(:,1) == t2);
        word_targ3_idx = -1;
        if numSuccessTargets == 2 || numSuccessTargets == 3
            word_targ3_idx = find(words(:,1) == t3);
        end
    else
        word_targ2_idx = -1;
        word_targ3_idx = -1;
    end

    % Determine previous end of trial index
    word_prev_end_idx = word_start_idx - 1;
    if word_prev_end_idx == 0
       word_prev_end_idx = -1;
    end
    % Return vector of word indices
    word_idx = [word_prev_end_idx word_start_idx word_targ1_idx word_targ2_idx word_targ3_idx word_end_idx];
end

% Currently there is something wrong with the extraction of words for the second '33' word trial. I dont know what the fuck is wrong with it but it must be fixed. Have fun. My guess is it has
% somethign to do with the numSuccessTargets. But could be further up. Idk. Almost there.