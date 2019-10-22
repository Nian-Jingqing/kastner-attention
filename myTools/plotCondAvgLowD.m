function plotCondAvgLowD(projected, sub_title)
% define legend
cond{1} = 'TL-V';
cond{2} = 'BL-V';
cond{3} = 'TR-V';
cond{4} = 'BR-V';
cond{5} = 'TL-V';
cond{6} = 'BL-V';
cond{7} = 'TR-V';
cond{8} = 'BR-V';
cond{9} = 'TL-H';
cond{10} = 'BL-H';
cond{11} = 'TR-H';
cond{12} = 'BR-H';
cond{13} = 'TL-H';
cond{14} = 'BL-H';
cond{15} = 'TR-H';
cond{16} = 'BR-H';

conds = [1, 9, 2, 10, 3, 11, 4, 12]; % hard code to sort conditions based on similar cue locations
cmap = lines();

% start plotting for the subplot
i = 1;
a = 1;
b = 1;
%axisLim = [-50 100 -20 40 -10 20];
%axisLimScale = [0.4, 0.8, 0.05, 0.8, 0.02, 0.02, 0.02];
for icond = conds
    lowD_this_cond = squeeze(projected(:,:, icond));
    h(i) = plot3(lowD_this_cond(:, 1), lowD_this_cond(:, 2), lowD_this_cond(:, 3), 'Color', (1-0.15*b)*cmap(a, :), 'LineWidth', 2, 'DisplayName', cond{icond});
    hold on
    plot3(lowD_this_cond(1, 1), lowD_this_cond(1, 2), lowD_this_cond(1, 3), 'o', 'MarkerFaceColor', (1-0.15*b)*cmap(a,:), 'MarkerEdgeColor', (1-0.15*b)*cmap(a,:))
    i = i+1;
    a = a+b-1;
    b = 1+rem(b,2);
end
%axis(gca, axisLim*axisLimScale(nday))
set(gca, 'view', [-29.1000, -6.0000]);
legend(h(1:8));
legend('Location', 'best')
title(sub_title);
