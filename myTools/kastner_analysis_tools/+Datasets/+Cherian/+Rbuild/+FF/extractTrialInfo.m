function [ trialStruct ] = extractTrialInfo( samp_indices, spikeInfo, targStruct, sampleRate, mad )
    trialStruct = struct();
    j = 1;
    for i=1:size(samp_indices,2)
        % get trial indices
        trial_indices = samp_indices{ i };
        % set intertrial start
        intertrial_start_idx = trial_indices( 1 );
        % set trial start
        trial_start_idx = trial_indices( 2 );
        % set intertrial end based off trial start
        intertrial_end_idx = trial_start_idx - 1;

        %intertrial_x = x(:,intertrial_start_idx:intertrial_end_idx);
        %intertrial_target = target(:,intertrial_start_idx:intertrial_end_idx);
        %intertrial_emg = emg(:,intertrial_start_idx:intertrial_end_idx);
        %intertrial_spikes = spikes(:,intertrial_start_idx:intertrial_end_idx);
        trialStruct( j ).expName = strrep( mad.header.matfilename, '.mat', '' );
        trialStruct( j ).trialNum = i;
        trialStruct( j ).trialType = Datasets.Cherian.Rbuild.Utils.TrialType.InterTrial;
        trialStruct( j ).numSuccessTargets = nan;
        trialStruct( j ).tIndices = [intertrial_start_idx nan nan nan intertrial_end_idx];
        trialStruct( j ).abstSidx = intertrial_start_idx;
        trialStruct( j ).abst1idx = nan;
        trialStruct( j ).abst2idx = nan;
        trialStruct( j ).abst3idx = nan;
        trialStruct( j ).abstEidx = intertrial_end_idx;
        trialStruct( j ).reltSidx = intertrial_start_idx - intertrial_start_idx + 1;
        trialStruct( j ).relt1idx = nan;
        trialStruct( j ).relt2idx = nan;
        trialStruct( j ).relt3idx = nan;
        trialStruct( j ).reltEidx = intertrial_end_idx - intertrial_start_idx + 1;

        trialStruct( j ).targPos = [ nan nan; nan nan; nan nan; ];
        
        %trialStruct( j ).x = intertrial_x;
        %trialStruct( j ).target = intertrial_target;
        %trialStruct( j ).emg = intertrial_emg;
        %trialStruct( j ).spikes = intertrial_spikes;
        trialStruct( j ).spikeInfo = spikeInfo;
        trialStruct( j ).sampleRate = sampleRate;
        j = j + 1;
        targ1_on_idx = trial_indices(3);
        targ2_on_idx = trial_indices(4);
        targ3_on_idx = trial_indices(5);
        trial_end_idx = trial_indices(6);
        %trial_x = x(:,trial_start_idx:trial_end_idx);
        %trial_target = target(:,trial_start_idx:trial_end_idx);
        %trial_emg = emg(:,trial_start_idx:trial_end_idx);
        %trial_spikes = spikes(:,trial_start_idx:trial_end_idx);
        trial_numSuccessTargets = targStruct(i).numSuccessfulTargets;

        trialStruct( j ).expName = strrep( mad.header.matfilename, '.mat', '' );        
        trialStruct( j ).trialNum = i;
        if trial_numSuccessTargets == 3
            trialStruct( j ).trialType = Datasets.Cherian.Rbuild.Utils.TrialType.Success;
        else
            trialStruct( j ).trialType = Datasets.Cherian.Rbuild.Utils.TrialType.Failure;
        end
        trialStruct( j ).numSuccessTargets = trial_numSuccessTargets;
        trialStruct( j ).tIndices = [trial_start_idx targ1_on_idx targ2_on_idx targ3_on_idx trial_end_idx];
        trialStruct( j ).abstSidx = trial_start_idx;
        trialStruct( j ).abst1idx = targ1_on_idx;
        trialStruct( j ).abst2idx = targ2_on_idx;
        trialStruct( j ).abst3idx = targ3_on_idx;
        trialStruct( j ).abstEidx = trial_end_idx;
        trialStruct( j ).reltSidx = trial_start_idx - trial_start_idx + 1;
        trialStruct( j ).relt1idx = targ1_on_idx - trial_start_idx + 1;
        trialStruct( j ).relt2idx = targ2_on_idx - trial_start_idx + 1;
        trialStruct( j ).relt3idx = targ3_on_idx - trial_start_idx + 1;
        trialStruct( j ).reltEidx = trial_end_idx - trial_start_idx + 1;
        
        trialStruct( j ).targPos = [ targStruct( i ).target{ 1 }{ 2 }; ...
                                     targStruct( i ).target{ 2 }{ 2 }; ...
                                     targStruct( i ).target{ 3 }{ 2 }; ];

        %trialStruct( j ).x = trial_x;
        %trialStruct( j ).target = trial_target;
        %trialStruct( j ).emg = trial_emg;
        %trialStruct( j ).spikes = trial_spikes;
        trialStruct( j ).spikeInfo = spikeInfo;
        trialStruct( j ).sampleRate = sampleRate;
        j = j + 1;
    end
end
