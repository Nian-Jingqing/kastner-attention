function [neuronStruct] = prepForWindowedCoherence(slidingStart, slidingEnd, slidingStep, window_size, spiking_struct, lfp_struct, alignInds, condType, cueLoc, rfLoc)
    nTrials_total = length(spiking_struct);
    nNeurons = size(spiking_struct(1).spikes,1);
    total_windows = (slidingEnd - slidingStart)/slidingStep + 1;
    for i = 1:total_windows
         window_center = slidingStart + (i - 1)*50;
         window_start = window_center - window_size/2 + 1; % to make sure that the window size is 300ms
         window_end = window_center + window_size/2;
         spiking_tensor = zeros(nTrials_total, nNeurons, window_size ); % nTrial x nNeurons x nTimes
         lfp_tensor = zeros(nTrials_total, nNeurons, window_size );
         for itrial = 1: nTrials_total
             spiking_tensor(itrial,:,:) = spiking_struct(itrial).spikes(:, (alignInds(i) + window_start):(alignInds(i) + window_end));
             lfp_tensor(itrial,:,:) = lfp_struct(itrial).lfp_filtered(:, (alignInds(i) + window_start):(alignInds(i) + window_end));
         end
         spiking_tensor = permute(spiking_tensor, [2, 1, 3]); % nNeurons x nTrials x nTimes
         lfp_tensor = permute(lfp_tensor, [2, 1, 3]);
         if strcmp(condType, 'InRF')
             selectedLocPerNeuron = rfLoc(:, 1);
         elseif strcmp(condType, 'OppositeRF')
             selectedLocPerNeuron = rfLoc(:, 2);
         end
         for n = 1:nNeurons
             isTrialThisCond = cueLoc == selectedLocPerNeuron(n);
             spiking_thisWindow = squeeze(spiking_tensor(n, isTrialThisCond, :));
             lfp_thisWindow = squeeze(lfp_tensor(n, isTrialThisCond, :));
             mean_lfp = mean(lfp_thisWindow, 2);
             centered_lfp = bsxfun(@minus, lfp_thisWindow, mean_lfp);
             centered_norm_lfp = centered_lfp./std(centered_lfp, 0, 2);
             neuronStruct(n).window(i).spiking = spiking_thisWindow;
             neuronStruct(n).window(i).lfp = centered_norm_lfp;
             neuronStruct(n).window(i).windowCenterTime = window_center;
         end
     end
end
