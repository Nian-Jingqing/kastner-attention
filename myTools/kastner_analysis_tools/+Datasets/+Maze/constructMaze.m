function constructMaze(trial)
p = trial.PARAMS;

bgcolor = zeros(1,3) + 0.92;
% make the background slightly dark
pos = [p.frameLeft p.frameBottom p.frameRight-p.frameLeft p.frameTop-p.frameBottom];
rbg = rectangle('Position', pos);
rbg.FaceColor = bgcolor;
rbg.EdgeColor = 'none';

% top frame
w=p.frameWidth;
pos = [p.frameLeft-w p.frameTop p.frameRight-p.frameLeft+2*w w];
rb(1) = rectangle('Position', pos);

% bottom frame
pos = [p.frameLeft-w p.frameBottom-w p.frameRight-p.frameLeft+2*w w];
rb(2) = rectangle('Position', pos);

% left frame
pos = [p.frameLeft-w p.frameBottom-w w p.frameTop-p.frameBottom+2*w];
rb(3) = rectangle('Position', pos);

% right frame
pos = [p.frameRight p.frameBottom-w w p.frameTop-p.frameBottom+2*w];
rb(4) = rectangle('Position', pos);

for nframe = 1:numel(rb)
    rb(nframe).FaceColor = zeros(1,3) + 0.3;
    rb(nframe).EdgeColor = 'none';
end


% plot some targets
for nf = trial.whichFly % 1:trial.numFlies
    w = p.flySize; % half width of the targets
    fx = p.flyX(nf)-w;
    fy = p.flyY(nf)-w;
    rt(nf) = rectangle('Position', [fx fy w*2 w*2]);
    rt(nf).FaceColor = [0.5 0 0];
    rt(nf).EdgeColor = 'none';
end

barrierColor = zeros(1,3) + 0.3;

% plot the barriers
for nb = 1:trial.numBarriers
    b = trial.BARRIER(nb);
    xstart = b.X - b.halfWidth;
    ystart = b.Y - b.halfHeight;
    pos = [xstart ystart b.halfWidth*2 b.halfHeight*2];
    bar(nb) = rectangle('Position', pos);
    bar(nb).FaceColor = barrierColor;
    bar(nb).EdgeColor = 'none';
end

% plot the start point
w = p.fixWindow;
pos = [p.fixX - w p.fixY-2 2*w 2*w];
rstart = rectangle('Position', pos);
rstart.FaceColor = zeros(1,3) + 1;
rstart.EdgeColor = zeros(1,3);

