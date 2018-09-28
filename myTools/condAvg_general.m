function condAvg = condAvg_general( data, fieldToCondAvg, condField )
% Summary of this function goes here
%   Detailed explanation goes here
% this function compute condition average and concatenate the condition-avg
% data together (nNeurons x (nConds x nTimes));

%% get total number of conditions
condVector = unique([data.( condField )]);
nCond = length(condVector);

%%  loop over all condition to compute condition avg for each condition
condStruct(nCond).avgInField = 1;

for cond = 1: nCond
    trialIndices = [ data.(condField) ] == condVector(cond);
    fieldTensor = permute(cat( 3, data(trialIndices).( fieldToCondAvg ) ), [ 3 2 1 ]);
    condStruct(cond).avgInField = squeeze(mean(fieldTensor, 1))';
end

%% concatenate the conditions together
condAvg = [condStruct.avgInField];

end