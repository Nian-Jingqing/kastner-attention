figure
for i = 1:16
    subplot(4,4,i)
    %plot(norm_rates(i,:))
    plot(tmp2(i,:))
    hold on
    plot(tmp(i,:))
    title(int2str(i))
end


%%

figure
i = 13;
plot(tmp2(i, :), 'LineWidth', 3)
min_s = min(tmp2(i,:));
max_s = max(tmp2(i,:));
height = max_s-min_s;
hold on

%
ns = find(tmp(i, :) > 0);
for j = 1:length(ns)
    y = [min_s, min_s+tmp(i, ns(j))*0.08*height];
    x = [ns(j), ns(j)];
    plot(x, y, 'Color', 'k', 'LineWidth', 6);
    hold on
end

ax=gca;
ax.XTick = [1 60 81];ax.XTickLabel={'-600','Target','+200'};
ax.YTick = [40-min_s];ax.YTickLabel={int2str(40-min_s)};
axis tight
xlabel('Time (ms)')
ylabel('Firing Rate (Spikes/s)')
title('Channel 13')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [234 369 1183 549])


%%
figure
imagesc(norm_rates);
ax=gca;
ax.XTick = [1 60 81];ax.XTickLabel={'-600','Target','+200'};
ax.YTick = [40-min_s];ax.YTickLabel={int2str(40-min_s)};
xlabel('Time (ms)')
ylabel('Channels')
title('Trial 150')
set(gca, 'FontSize', 20)
set(gcf, 'Position', [100 230 1136 374])