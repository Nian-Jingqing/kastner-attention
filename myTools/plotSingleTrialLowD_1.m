function plotSingleTrialLowD_1(alf, UE, nday, pca_proj_mat, rateOrSpike, timePoints, sub_title, all_data_means)
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
%for nCue = 1:numel(conditions.cueOn)
for nCue = [2, 4]
    %for nBar = 1:numel(conditions.barOn)
    for nBar = 1
        trialIndicesThisCond = find((UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue));
        for itr = 1:numel(trialIndicesThisCond)
            ntr = trialIndicesThisCond(itr);
            data_this_trial = alf{nday}(ntr).(rateOrSpike)(:, alf{nday}(ntr).cueOnset + timePoints);
            data_this_trial_centered = bsxfun(@minus, data_this_trial, mean(data_this_trial, 2));
            %data_this_trial_centered = bsxfun(@minus, data_this_trial, all_data_means);
            lowD_this_trial = pca_proj_mat' * data_this_trial_centered;
            
            if itr == 1
                %h(i) = plot3(lowD_this_trial(1, :), lowD_this_trial(2, :), lowD_this_trial(3, :), 'Color', (1-0.15*b)*cmap(a, :), 'LineWidth', 1, 'DisplayName', cond{i});
                h(i) = plot3(lowD_this_trial(1, :), lowD_this_trial(2, :), lowD_this_trial(3, :), 'Color', (1-0.15*b)*cmap(nCue, :), 'LineWidth', 1, 'DisplayName', cond{2*nCue - 1});
            else
                %plot3(lowD_this_trial(1, :), lowD_this_trial(2, :), lowD_this_trial(3, :), 'Color', (1-0.15*b)*cmap(a, :), 'LineWidth', 1);
                plot3(lowD_this_trial(1, :), lowD_this_trial(2, :), lowD_this_trial(3, :), 'Color', (1-0.15*b)*cmap(nCue, :), 'LineWidth', 1);
            end
            
            hold on
            %plot3(lowD_this_trial(1, 1), lowD_this_trial(2, 1), lowD_this_trial(3, 1), 'o', 'MarkerFaceColor', (1-0.15*b)*cmap(a,:), 'MarkerEdgeColor', (1-0.15*b)*cmap(a,:))
            plot3(lowD_this_trial(1, 1), lowD_this_trial(2, 1), lowD_this_trial(3, 1), 'o', 'MarkerFaceColor', (1-0.15*b)*cmap(nCue,:), 'MarkerEdgeColor', (1-0.15*b)*cmap(nCue,:))
        end
        i = i+1;
        a = a+b-1;
        b = 1+rem(b,2);
    end
end
set(gca, 'view', [-40.7000, -22.8000]);
legend(h(1:2));
%legend(h(1:8));
legend('Location', 'best')
title(sub_title);