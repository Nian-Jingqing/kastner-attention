function [spikeBandData] = lowPassNotchBandPass(dataBuffer_normd)
% first apply low pass
Fs = 40000;
% apply low pass first
[b,a] = butter(4, [5000] / (Fs / 2), 'low');
dataBuffer_normd = filtfilt(b, a, dataBuffer_normd');
dataBuffer_normd = dataBuffer_normd';

% then perform pwelch and find peaks
Pxx_postLow = [];
for nn = 1:size(dataBuffer_normd, 1)
    [ Pxx_postLow(:,nn), w ] = pwelch( dataBuffer_normd( nn, : ), 40000*4 );
end

x = w * Fs/(2*pi);
locs = {};
pks = {};
for nn = 1:size(dataBuffer_normd, 1)
    [pks_tmp, locs_tmp] = findpeaks(10 * log10(Pxx_postLow(:, nn)'), x,  'MinPeakProminence', 5 );
    %pksInBand = find(locs_tmp <= 5000 & locs_tmp >= 300);
    locs{nn} = locs_tmp;
    pks{nn} = pks_tmp;
end
dataFiltered = dataBuffer_normd;

% then perform notch filter
parfor nn = 1:size(dataFiltered, 1)

    % % dummy filter
    %hcas = dfilt.df1( 1, 1);
    for pk = 1:numel(locs{nn})
        freqToNotch = locs{nn}(pk); 
        W0 = freqToNotch/(Fs/2); BW = 4/(Fs/2); %BW = W0/35;
        [b,a] = iirnotch(W0, BW);

        % cascade filters
        %hnew = dfilt.df1( b, a);
        %hcas = dfilt.cascade( hcas, hnew );

        % old filtering code
        dataFiltered(nn, :) = filter(b,a,dataFiltered(nn,:));
    end

    % % cascade filters
    %dataFiltered(nn, :) = hcas.filter( dataFiltered(nn,:) );
    
    %disp(['finished ' int2str(pk)])
end

% then perform bandpass filtering
filtLowCutoff = 300;
filtHighCutoff = 5000;
[b,a] = butter(4, [filtLowCutoff filtHighCutoff] / (Fs / 2), 'bandpass');
spikeBandData = filtfilt(b, a, dataFiltered');

spikeBandData = spikeBandData';

