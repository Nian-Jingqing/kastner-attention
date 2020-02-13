function [cond_avg] = getCondAvgAllDays_factors(dataset, UE, alf, timePoints, nTimes, nConds, conditions, binsize_rescaled)
% cond_avg is the condition avg matrix (nChannel x nConditions*nTimes).
% There are 2 bar types and 8 cue types, so combined 16 conditions at cueOnset
% 16 is hard coded -- need to change
num_factors = size(alf{1}(1).factors, 1);
times_base = 1:nTimes;
cond_avg = zeros(num_factors, nConds*nTimes);

for nBar = 1:numel(conditions.barOn)
    for nCue = 1:numel(conditions.cueOn)
        cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
        trial_num = 0;
        this_cond_sum = zeros(num_factors, nTimes);
        for nday = 1:numel(dataset)
            trialIndices = find((UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue));
            for itr = 1:numel(trialIndices)
                ntr = trialIndices(itr);
                this_cond_sum = this_cond_sum + alf{nday}(ntr).factors(:, alf{nday}(ntr).cueOnset + timePoints);
            end
            trial_num = trial_num + numel(trialIndices);
        end
        this_cond_inds = ((nBar - 1)*8 + nCue - 1) * nTimes + times_base;
        cond_avg(:, this_cond_inds) = this_cond_sum / trial_num;
    end
end
