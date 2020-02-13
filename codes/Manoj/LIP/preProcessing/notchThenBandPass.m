function [spikeBandData] = notchThenBandPass(dataBuffer_normd)
% first apply low pass
Fs = 40000;
% then perform pwelch and find peaks
Pxx_normd = [];
for nn = 1:size(dataBuffer_normd, 1)
    [ Pxx_normd(:,nn), w ] = pwelch( dataBuffer_normd( nn, : ), 40000*4 );
end

x = w * Fs/(2*pi);
locs = {};
pks = {};
for nn = 1:size(dataBuffer_normd, 1)
    [pks_tmp, locs_tmp] = findpeaks(10 * log10(Pxx_normd(:, nn)'), x,  'MinPeakProminence', 5 );
    pksInBand = find(locs_tmp <= 5010 & locs_tmp >= 290);
    locs{nn} = locs_tmp(pksInBand);
    pks{nn} = pks_tmp(pksInBand);
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

