function out = calcXcorr(seq, timeWindow, neuralField, lfpField, neuralResample, numRandom, whichChannels)
out = [];

numChannels = size(seq(1).(neuralField), 1);

if ~exist('neuralResample','var') || isempty(neuralResample)
    neuralResample = 1;
end

if ~exist('numRandom','var')
    numRandom = 2;
end

if ~exist('whichChannels', 'var')
    whichChannels = 1:numChannels;
end


maxlag = floor( diff ( timeWindow ) /2 );

% set the rand number generator
rng('default')

for nch = 1:numChannels
    if ~intersect(whichChannels, nch)
        continue;
    end

    clear xvals xinds randxvals whichTrial
    trials_used = 0;
    % initially allocate a large matrix, then trim unused trials (because of NANs)
    xvals = zeros(numel(seq), maxlag *2+1);
    for ntr = 1:numel(seq)
        y = seq(ntr).(neuralField)(nch, :);
        y = processy(y, neuralResample, timeWindow);

        lfp = seq(ntr).LFP(nch, timeWindow(1):timeWindow(2));
        lfp = lfp - mean(lfp);
        if any(isnan(lfp(:)))
            continue;
        end
        trials_used = trials_used + 1;

        [a,b] = xcorr( full(y), lfp, maxlag);

        xvals(trials_used, :) = a;
        xinds = b;
        whichTrial(trials_used) = ntr;
    end
    out(nch).xvals = xvals(1:trials_used, :);
    out(nch).xinds = xinds;
    out(nch).whichTrial = whichTrial;


    allRandxvals = zeros(numel(seq), maxlag *2+1);
    %allRandxvals = zeros(numRandom, maxlag*2 + 1);
    trials_used = 0;
    randxvals = zeros(1, maxlag*2 + 1);
    for ntr = 1:numel(seq)
        y = seq(ntr).(neuralField)(nch, :);
        y = processy(y, neuralResample, timeWindow);

        randtrial = floor(rand() * numel(seq))+1;
        lfp = seq(randtrial).LFP(nch, timeWindow(1):timeWindow(2));
        while(any(isnan(lfp(:))) || ntr == randtrial)
            randtrial = floor(rand() * numel(seq))+1;
            lfp = seq(randtrial).LFP(nch, timeWindow(1):timeWindow(2));
        end
        lfp = lfp - mean(lfp);
        if any(isnan(lfp(:)))
            continue;
        end
        trials_used = trials_used + 1;
        [randa,randb] = xcorr( full(y),lfp, maxlag);

        allRandxvals(trials_used, :) = randa;
        whichTrial(trials_used) = ntr;
    end


    out(nch).randxvals = allRandxvals(1:trials_used, :);
end

function yout = processy(y, neuralResample, timeWindow)
y = y-mean(y);
% if needed, resample y
if neuralResample ~= 1
    y = resample(y', neuralResample, 1)';
end
    yout = y( timeWindow(1):timeWindow(2) );