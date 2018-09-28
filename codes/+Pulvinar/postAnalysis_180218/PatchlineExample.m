t = 0:pi/64:4*pi;
p(1) = Plot.patchline(t,sin(t),'edgecolor','b','linewidth',2,'edgealpha',0.5);
p(2) = Plot.patchline(t,cos(t),'edgecolor','r','linewidth',2,'edgealpha',0.5);
l = legend('sine(t)','cosine(t)');
tmp = sort(findobj(l,'type','patch'));
for ii = 1:numel(tmp)
    set(tmp(ii),'facecolor',get(p(ii),'edgecolor'),'facealpha',get(p(ii),'edgealpha'),'edgecolor','none')
end