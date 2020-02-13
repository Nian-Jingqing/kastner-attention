function plotSingleTrialPC(cueID, pcID, alf, UE, nday, proj_matrix_this_day, rateOrSpike, timePoints, sub_title, proj_bias_this_day, binsize_rescaled)
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
hold all
for nCue = cueID
    %for nBar = 1:numel(conditions.barOn)
    for nBar = 1
        trialIndicesThisCond = find((UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue));
        for itr = 1:numel(trialIndicesThisCond)
            ntr = trialIndicesThisCond(itr);
            if strcmp(rateOrSpike, 'spikes')
                data_this_trial = (1000/binsize_rescaled)*alf{nday}(ntr).(rateOrSpike)(:, alf{nday}(ntr).cueOnset + timePoints);
            else
                
                data_this_trial = alf{nday}(ntr).(rateOrSpike)(:, alf{nday}(ntr).cueOnset + timePoints);
            end
            %data_this_trial_centered = bsxfun(@minus, data_this_trial, mean(data_this_trial, 2)); % not correct
            data_this_trial_centered = bsxfun(@minus, data_this_trial, proj_bias_this_day);
            lowD_this_trial = proj_matrix_this_day' * data_this_trial_centered; % this should be nChannels x nTimes
            %sigma = 3;
            %lowD_this_trial = smooth2Dmatrix(lowD_this_trial, sigma);
            if itr == 1
                h(i) = plot(lowD_this_trial(pcID, :), 'Color', 0.85*cmap(cueID, :), 'LineWidth', 1, 'DisplayName', sub_title);
            else
                plot(lowD_this_trial(pcID, :), 'Color', 0.85*cmap(cueID, :), 'LineWidth', 1);
            end            
            hold on
        end
    end
end
%legend(h(1));
