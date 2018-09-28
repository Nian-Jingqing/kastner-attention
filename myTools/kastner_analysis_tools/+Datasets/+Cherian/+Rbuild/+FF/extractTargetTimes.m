function [target_times] = extractTargetTimes(targStruct,i)
    t1 = targStruct(i).target{1}{1};
    t2 = targStruct(i).target{2}{1};
    t3 = targStruct(i).target{3}{1};
    target_times = [t1 t2 t3];
    
end
