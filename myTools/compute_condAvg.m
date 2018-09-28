function [ condAvg, condAvg_up, condAvg_low ] = compute_condAvg( RstructCopy, fieldToCondAvg, condType )
% Summary of this function goes here
%   Detailed explanation goes here
%% get the trials for the selected alignment type
condIx = [RstructCopy.r.conditionId]; 
% holdIx = [RstructCopy.r.isHoldTrial];
realAlign = RstructCopy.r(condIx == condType);
nTrials = length(realAlign);


%% calculate condition-avg

fieldTensor = permute(cat( 3, realAlign.( fieldToCondAvg ) ), [ 3 2 1 ]);
% extract out the trials for the condition and re-arrange the matrix
% to nTrial x nTimes x nThings
condAvg = squeeze(mean(fieldTensor, 1))';
std_field = squeeze(std(fieldTensor, 1))';
stderr_field = std_field./(sqrt(nTrials));
condAvg_up = condAvg + stderr_field;
condAvg_low = condAvg - stderr_field;


end