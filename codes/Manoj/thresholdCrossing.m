%filename = '/mnt/scratch/feng/Remy_02262019_PUL_Raw.pl2';
%outfile = '/mnt/scratch/cpandar/Remy_02262019_PUL_spikeband.mat';
%
%tic;
%broadband2streamMinMax( filename, outfile )
%toc;

%%
spikeBandFile = '/snel/share/share/data/kastner/Manoj/PUL/spikeBand/Remy_02182019_PUL_spikeband.mat';
bb = load(spikeBandFile);
bb = bb.spikeband;

%% get var and remove NaN, also remove NaN for minSpikeBand
for ich = 1:size(bb.minSpikeBand,2)
    chVar{ ich } = bb.meanSquared( bb.meanSquaredChannel == ich );
    whereNan = find( isnan( chVar{ ich } ) );
    chVar{ ich } = chVar{ ich }( 1 : (whereNan(1) - 1) );
    whereNan_msb = find( isnan( bb.minSpikeBand( : , ich ) ) );
    chMsb{ ich } = bb.minSpikeBand(1:(whereNan_msb(1) - 1), ich);
end

%% compute changing std
chStdVec = zeros(1,32);
for ich = 1:numel( chVar )
    % % use original variance (changing value across a session)
    chStd{ich} = sqrt(chVar{ich});
    chStd{ich} = repelem(chStd{ich}, 32);
    if numel(chStd{ich}) < numel(chMsb{ich})
        num_diff = numel(chMsb{ich}) - numel(chStd{ich});
        elemToAdd = chStd{ich}(end,1)*ones(num_diff, 1);
        chStd{ich} = [chStd{ich};elemToAdd];
    else
        chStd{ich} = chStd{ich}(1:length(chMsb{ich}));
    end 
end

%% compute constant std
chStdVec = zeros(1,32);
for ich = 1:numel( chVar )
    %use mean of variance across entire session for each channel
    chStd_cons{ich} = sqrt( mean( chVar{ich} ) );
    chStdVec(ich) = chStd_cons{ich};
end

%% remove common average
bb.minSpikeBand_rmCA = bb.minSpikeBand - mean(bb.minSpikeBand,2);

%% get spikes
spikes = sparse(size(chMsb{1}, 1), numel(chMsb));
for ich = 1:size(bb.minSpikeBand,2)
    chMean = mean(chMsb{ich});
    chThres = chMean - 2*chStd_cons{ich};
    leftValue = chMsb{ich} - chThres;
    spikes(leftValue <= 0, ich) = 1;
end

% remove spikes if show on more than 25 TCs (~80%)
allSpikesMS = sum(spikes,2);
spikes((allSpikesMS > 25), :) = 0;

% get rid of MU 1 - 4
spikes = spikes(:, 5:end);
stream.spikes = sparse(spikes);
%%
verifyThresholdCrossings
