function plotSingleTrialLowD(alf, UE, nday, proj_matrix_this_day, rateOrSpike, timePoints, sub_title)
% define legend
cond{1} = 'TL-V';
cond{2} = 'TL-H';
cond{3} = 'BL-V';
cond{4} = 'BL-H';
cond{5} = 'TR-V';
cond{6} = 'TR-H';
cond{7} = 'BR-V';
cond{8} = 'BR-H';

%conds = [1, 9, 2, 10, 3, 11, 4, 12]; % hard code to sort conditions based on similar cue locations
cmap = lines();

conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR'};
trialIndicesThisCond = [];
i = 1;
a = 1;
b = 1;
%keyboard
for nCue = 1:numel(conditions.cueOn)
    for nBar = 1:numel(conditions.barOn)
        trialIndicesThisCond = find((UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue));
        for itr = 1:numel(trialIndicesThisCond)
            ntr = trialIndicesThisCond(itr);
            data_this_trial = alf{nday}(ntr).(rateOrSpike)(:, alf{nday}(ntr).cueOnset + timePoints);
            data_this_trial_centered = bsxfun(@minus, data_this_trial, mean(data_this_trial, 2));
            lowD_this_trial = proj_matrix_this_day' * data_this_trial_centered; % this should be nChannels x nTimes
            
            if itr == 1
                h(i) = plot3(lowD_this_trial(1, :), lowD_this_trial(2, :), lowD_this_trial(3, :), 'Color', (1-0.15*b)*cmap(a, :), 'LineWidth', 1, 'DisplayName', cond{i});
            else
                plot3(lowD_this_trial(1, :), lowD_this_trial(2, :), lowD_this_trial(3, :), 'Color', (1-0.15*b)*cmap(a, :), 'LineWidth', 1);
            end
            
            hold on
            plot3(lowD_this_trial(1, 1), lowD_this_trial(2, 1), lowD_this_trial(3, 1), 'o', 'MarkerFaceColor', (1-0.15*b)*cmap(a,:), 'MarkerEdgeColor', (1-0.15*b)*cmap(a,:))
        end
        i = i+1;
        a = a+b-1;
        b = 1+rem(b,2);
    end
end
set(gca, 'view', [-40.7000, -22.8000]);
legend(h(1:8));
legend('Location', 'best')
title(sub_title);