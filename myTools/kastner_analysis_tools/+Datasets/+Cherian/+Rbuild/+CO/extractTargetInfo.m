function [ targetStruct ] = extractTargetInfo( targ_table, mA )
    cell_targets = table2cell(targ_table);
    analogStartTime = mA(1).info.starttime;
    targetStruct = struct();
    for i=1:size(cell_targets,1)
        % Find target Onset Time
        
        targOnTime = cell_targets{i,1}(1) - analogStartTime;
        goCueTime = cell_targets{i,2} - analogStartTime;
        
        Xpos = cell_targets{i,1}(2); 
        Ypos = cell_targets{i,1}(3); 
            
        targPosition = [Xpos Ypos];
        
        targetStruct(i).pos = targPosition;
        targetStruct(i).targOnTime = targOnTime;
        targetStruct(i).goCueTime = goCueTime;
    end
end
