function f = blankFigure( fnum )

% use an existing figure, or not?
if exist('fnum', 'var' )
    f = figure(fnum);
    clf;
else
    f = figure();
end

set(gcf,'color',[1 1 1]);
h = subplot('Position', [0 0 1 1]);
set(h, 'visible', 'off');