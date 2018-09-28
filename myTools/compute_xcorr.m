function [xcorr_struct ] = compute_xcorr(nNeurons, real_allConds_dataStruct, lfads_allConds_dataStruct, real_chopDuration, lfads_chopDuration, maxLag, par, rfLoc, cueLoc, whichType )
% the difference between real and lfads is that lfads data is rebinned and real has the raw sampling rate
    xcorr_struct = struct;
    for n = 1:nNeurons
        if strcmp(whichType, 'inRF')
            real_dataStruct = real_allConds_dataStruct(cueLoc == rfLoc(n, 1));
            lfads_dataStruct = lfads_allConds_dataStruct(cueLoc == rfLoc(n, 1));
        elseif strcmp(whichType, 'oppositeRF')
            real_dataStruct = real_allConds_dataStruct(cueLoc == rfLoc(n, 2));
            lfads_dataStruct = lfads_allConds_dataStruct(cueLoc == rfLoc(n, 2));
            end
        nTrials = length(real_dataStruct);
        allTrialIndices = 1: nTrials;
        for itrial = 1: nTrials

            
            % spiking = real_dataStruct(itrial).spiking(n,1:real_chopDuration);
            % allOtherTrialIndices = allTrialIndices(allTrialIndices~=itrial);
            % oneRandomTrial = randsample(allOtherTrialIndices,1);
            % lfp_spiking = real_dataStruct(itrial).lfp(n,1:real_chopDuration);
            % lfp_spiking_shuffled = real_dataStruct(oneRandomTrial).lfp(n,1:real_chopDuration);
            % xcorr_struct(n).r_spikingCorrelation(itrial,:) = xcorr(spiking, lfp_spiking, maxLag);
            % xcorr_struct(n).r_spikingCorrelation_shuffled(itrial,:) = xcorr(spiking, lfp_spiking_shuffled, maxLag);

            spiking = lfads_dataStruct(itrial).lfadsRate(n,1:lfads_chopDuration);
            allOtherTrialIndices = allTrialIndices(allTrialIndices~=itrial);
            oneRandomTrial = randsample(allOtherTrialIndices,1);
            lfp_spiking = lfads_dataStruct(itrial).lfp(n,1:lfads_chopDuration);
            lfp_spiking_shuffled = lfads_dataStruct(oneRandomTrial).lfp(n,1:lfads_chopDuration);
            xcorr_struct(n).r_spikingCorrelation(itrial,:) = xcorr(spiking, lfp_spiking, maxLag/par.spikeBinMs);
            xcorr_struct(n).r_spikingCorrelation_shuffled(itrial,:) = xcorr(spiking, lfp_spiking_shuffled, maxLag/par.spikeBinMs);

            lfadsRate = lfads_dataStruct(itrial).lfadsRate(n,1:lfads_chopDuration);
            lfp_lfads = lfads_dataStruct(itrial).lfp(n,1:lfads_chopDuration);
            lfp_lfads_shuffled = lfads_dataStruct(oneRandomTrial).lfp(n,1:lfads_chopDuration);
            xcorr_struct(n).r_lfadsCorrelation(itrial,:) = xcorr(lfadsRate, lfp_lfads, maxLag/par.spikeBinMs);
            xcorr_struct(n).r_lfadsCorrelation_shuffled(itrial,:) = xcorr(lfadsRate, lfp_lfads_shuffled, maxLag/par.spikeBinMs);
        end
    end
end
