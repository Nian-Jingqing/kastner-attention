function[samp_indices] = generateTrialSampleIndices(targStruct,mA,words,sampleRate)
    for i=1:size(targStruct,2)
        % Number of succesful targets reached during trial
        numSuccessTargets = targStruct(i).numSuccessfulTargets;
        % Number of targets that appeared during trial 
        numTargetsAppear = numSuccessTargets +1;
        % Extract times when targets appeared during trial
        target_times = Datasets.Cherian.Rbuild.FF.extractTargetTimes( targStruct, i );

        % Extract indices in words table of relevenat trial times
        word_indices = Datasets.Cherian.Rbuild.FF.extractWordIndices(words,target_times,numSuccessTargets);

        % Extract sample indices for relevant trial times
        tmp_samp_indices = Datasets.Cherian.Rbuild.FF.extractSampleIndicesFromWordIndices(word_indices,words,mA,sampleRate);

        % Push sample indices into a cell holding sample indices for all trials
        samp_indices{i} = tmp_samp_indices;
    end
end
