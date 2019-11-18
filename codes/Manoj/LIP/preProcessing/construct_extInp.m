function [extInp] = construct_extInp(UE, total_samples)
%% construct external inputs
UE.cueEnd = UE.cueOn + 150;
UE.targetEnd = UE.targetOn + 150;
UE.targetEnd(UE.targetEnd >= UE.trialEnd) = UE.trialEnd(UE.targetEnd >= UE.trialEnd);
% make sure the targetEnd is the same as the trialEnd for the early stopping trials

nCueTypes = numel(unique(UE.cueType));
nFixTypes = numel(unique(UE.fixType));
nBarTypes = numel(unique(UE.barType));
nTargetTypes = 4; 

fixDims = sparse(total_samples, 1 + nFixTypes);
barDims = sparse(total_samples, nBarTypes);
cueDims = sparse(total_samples, nCueTypes/2); % because these cueDims are only for Exo
targetDims = sparse(total_samples, nTargetTypes);
% construct fix square dimensions
% Total five dimensions. First dimension - plain color. The rest four dimensions - cue colors

for i = 1:numel(UE.fixOn)
    
    % construct fix square dimensions
    % Total five dimensions. First dimension - plain color. The rest four dimensions - cue colors
    if UE.fixType(i) == 1
        fixDims(UE.fixOn(i) + 35 : UE.trialEnd(i) + 35, 1) = 1; % first dimension is plain color
    elseif UE.fixType(i) == 2 % for color cue, turn down first dim during the cue
        fixDims(UE.fixOn(i) + 35 : UE.cueOn(i) + 35, 1) = 1;
        fixDims(UE.cueEnd(i) + 35 : UE.trialEnd(i) + 35, 1) = 1;
        fixDims(UE.cueOn(i) + 35 : UE.cueEnd(i) + 35, UE.cueType(i) - (nCueTypes/2 - 1)) = 1;
    end

    % construct bar dimensions

    barDims(UE.barOn(i) + 35 : UE.trialEnd(i) + 35, UE.barType(i)) = 1;

    % construct cue dimensions (only for exo cues, because endo are changes in fix))
    if UE.fixType(i) == 1
        cueDims(UE.cueOn(i) + 35 : UE.cueEnd(i) + 35, UE.cueType(i)) = 1;
    end

    % construct target dimensions
    targetPos = Manoj_getTargetPos(UE.cueType(i));
    if UE.isValidTarget(i)
        targetDims(UE.targetOn(i) + 35 : UE.targetEnd(i) + 35, targetPos) = 1;
    elseif (UE.isSameObjTarget(i) && UE.barType(i) == 1) || (UE.isDiffObjTarget(i) && UE.barType(i) == 2)
        if ismember(targetPos, [1, 3])
            shiftedTargetPos = targetPos + 1;
        else
            shiftedTargetPos = targetPos - 1;
        end
        targetDims(UE.targetOn(i) + 35 : UE.targetEnd(i) + 35, shiftedTargetPos) = 1;
    elseif (UE.isSameObjTarget(i) && UE.barType(i) == 2) || (UE.isDiffObjTarget(i) && UE.barType(i) == 1)
        if ismember(targetPos, [1, 2])
            shiftedTargetPos = targetPos + 2;
        else
            shiftedTargetPos = targetPos - 2;
        end        
        targetDims(UE.targetOn(i) + 35 : UE.targetEnd(i) + 35, shiftedTargetPos) = 1;
    end

end
%%
extInp = [fixDims, barDims, cueDims, targetDims];