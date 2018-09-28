function [ condAvgFiringRate ] = condAvg_onSelectedAlignType( r_realCopy, alignType, isHold, binSize )
% Summary of this function goes here
%   Detailed explanation goes here
%% get the trials for the selected alignment type
alignIx = [r_realCopy.r.alignType]; 
holdIx = [r_realCopy.r.isHoldTrial];
realAlign = r_realCopy.r(alignIx == alignType & holdIx == isHold);


%% initialize a struct array to store avg firing rate for each condition (cue loc)
nCond = length(unique([realAlign.cueLoc]));
% get condition number
condStruct(nCond).avgSmoothedSpikes = 1;
% pre-allocate the condStruct to avoid changing size in the for loop

for cond = 1: nCond
    trialIndices = [ realAlign.cueLoc ] == cond;
    spikeMatrixPerCond = permute(cat(3, realAlign(trialIndices).smoothed_cutOff), [ 3 2 1 ]);
    % extract the trials for the condition out and re-arrange the matrix
    % to nTrial x nTimes x nNeurons
    condStruct(cond).avgSmoothedSpikes = squeeze(mean(spikeMatrixPerCond, 1))';
    % calculate avg firing rate for each condition, nNeurons x nTimes
    
end

condAvgFiringRate = [condStruct.avgSmoothedSpikes];
condAvgFiringRate = condAvgFiringRate * (1000/binSize);
% %% get avg firing rate for all the condition and substract this avg firing rate for all the conditions
% 
% spikeTensorAllConds = permute(cat(3, condStruct.avgSmoothedSpikes), [ 3 2 1 ]);
% % get a spike tensor for all cond: nCond x nTimes x nNeurons
% avgFiringRateAllConds = squeeze(sum(spikeTensorAllConds, 1)/nCond)';
% % calculate the avg firing rate for all conditions, nNeurons x nTimes
% 
% for cond = 1:nCond
%     condStruct(cond).avgSmoothedSpikes_Norm = condStruct(cond).avgSmoothedSpikes - avgFiringRateAllConds;
%     % normalize the avg smoothed spikes for each condition
% end   
% condAvgFiringRate = [condStruct.avgSmoothedSpikes_Norm];
% condAvgFiringRate = condAvgFiringRate * (1000/binSize);
% % concatenate the normalized cond-avg firing rate together, nNeurons x (nTimes x nCond)

end

