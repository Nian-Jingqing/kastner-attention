function [ sample_idx ] = extractSampleIndicesFromWordIndices( word_idx, words, mA, sampleRate)
    word_prev_end_idx = word_idx(1);
    word_start_idx = word_idx(2);
    word_targ1_idx = word_idx(3);
    word_targ2_idx = word_idx(4);
    word_targ3_idx = word_idx(5);
    word_end_idx = word_idx(6);

    % Extract trial start index
    trial_start_idx = Datasets.Cherian.Rbuild.Utils.convertWordIdx2sampleIdx(word_start_idx,words,sampleRate);

    % Extract target 1 index
    targ1_on_idx = Datasets.Cherian.Rbuild.Utils.convertWordIdx2sampleIdx(word_targ1_idx,words,sampleRate);

    % If target 2 exists then find index
    if word_targ2_idx == -1
        targ2_on_idx = nan;
    else
        targ2_on_idx = Datasets.Cherian.Rbuild.Utils.convertWordIdx2sampleIdx(word_targ2_idx,words,sampleRate);
    end
    
    % If target 3 exists then find index    
    if word_targ3_idx == -1
        targ3_on_idx = nan;
    else
        targ3_on_idx = Datasets.Cherian.Rbuild.Utils.convertWordIdx2sampleIdx(word_targ3_idx,words,sampleRate);
    end

    % If end trial exists find index, if not set it as the end
    if word_end_idx == -1
        word_end_idx = size(words,1);
        trial_end_idx = Datasets.Cherian.Rbuild.Utils.convertWordIdx2sampleIdx(word_end_idxk,words,sampleRate);
    else
        trial_end_idx = Datasets.Cherian.Rbuild.Utils.convertWordIdx2sampleIdx(word_end_idx,words,sampleRate);
    end

    % If previous end of trial exists find index, if not set it as the first sample
    if word_prev_end_idx == -1
        prev_trial_end_idx = 1;
    else
        prev_trial_end_idx = Datasets.Cherian.Rbuild.Utils.convertWordIdx2sampleIdx(word_prev_end_idx,words,sampleRate)+1;
    end
    
    sample_idx = [prev_trial_end_idx trial_start_idx targ1_on_idx targ2_on_idx targ3_on_idx trial_end_idx];
end
