
%% low pass
dataBuffer_normd = normalize(dataBuffer, 'centered');
Fs = 40000;
% apply low pass first
[b,a] = butter(10, [5000] / (Fs / 2), 'low');
dataBuffer_normd = filtfilt(b, a, dataBuffer_normd');
dataBuffer_normd = dataBuffer_normd';
%%

%
Fs = 40000;
Pxx_normd = [];
for nn = 1:size(dataBuffer_normd, 1)
    [ Pxx_normd(:,nn), w ] = pwelch( dataBuffer_normd( nn, : ), 40000*4 );
end

x = w * Fs/(2*pi);
locs = {};
pks = {};
for nn = 1:size(dataBuffer_normd, 1)
    [pks_tmp, locs_tmp] = findpeaks(10 * log10(Pxx_normd(:, nn)'), x,  'MinPeakProminence', 5 );
    %pksInBand = find(locs_tmp <= 5000 & locs_tmp >= 300);
    locs{nn} = locs_tmp;
    pks{nn} = pks_tmp;
end

dataFiltered = dataBuffer_normd;



%%
Pxx_normd_2 = [];
for nn = 1:size(dataBuffer_normd, 1)
    [ Pxx_normd_2(:,nn), w ] = pwelch( dataBuffer_normd( nn, : ), 40000*4 );
end


%%
parfor nn = 1:size(dataBuffer_normd, 1)

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

%%
Pxx_normd_2 = [];
for nn = 1:size(dataBuffer_normd, 1)
    [ Pxx_normd_2(:,nn), w ] = pwelch( dataFiltered( nn, : ), 40000*4 );
end

x = w * Fs/(2*pi);

%% compare

clf
subplot(2,1,1)
plot(x, 10*log10(Pxx_normd(:,1)'))
%plot(x, mean(10*log10(Pxx_normd')))
xlabel( 'Freq (Hz)');
xlim([300 7000]);
%ylim([-20 60])
title('Pre SP')

subplot(2,1,2)
plot(x, 10*log10(Pxx_normd_2(:,1)'))
%plot(x, mean(10*log10(Pxx_normd_2')))
xlabel( 'Freq (Hz)');
xlim([300 7000]);
%ylim([-20 60])
title('Post SP')









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% comparing stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalize power
preSP_normd = normalize(dataBuffer, 'centered');
%clear dataBuffer

%% pwelch on normalized dataBuffer
Fs = 40000;

Pxx_preSP = [];
for nn = 1:32
    [ Pxx_preSP(:,nn), w ] = pwelch( preSP_normd( nn, : ), 40000*4 );
end

clear preSP_normd

%% plot
figure()
%ax(1) = subplot(2,1,1)
plot( w * Fs/(2*pi), mean(10 * log10( Pxx_preSP' )) );
xlabel( 'Freq ( Hz )' );
xlim([300 5000 ] );
title('Pre_signalProcessing')



%% load processed code
spikeBandBase = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/notchFilt_bandPass/tmp/';
%spikeBandBase = '/mnt/scratch/feng/LIP_spikeBand/';
spikeBandFileName = ['Remy_', dataset(day).date, '_LIP_spikeband.mat'];
spikeBandFile = fullfile(spikeBandBase, spikeBandFileName);
bb = load(spikeBandFile);
bb = bb.spikeband;


%%
NUM_SAMPLES_PLEXON_BLOCK = 65535; % PL2 data block size
MILLISECONDS_TO_GET = 120 * 1000;
% make that a multiple of the block size to simplify things
MILLISECONDS_TO_GET = ceil( MILLISECONDS_TO_GET / NUM_SAMPLES_PLEXON_BLOCK ) * NUM_SAMPLES_PLEXON_BLOCK;

% how many plexon samples per millisecond?
SAMPLES_PER_MS = Freq / 1000;
% make sure this is an int
assert( mod( SAMPLES_PER_MS, 1 ) == 0, 'not an integer number of samples' );

% Make data buffer this size
SAMPLES_IN_BUFFER = MILLISECONDS_TO_GET * SAMPLES_PER_MS;




%%
SAMPLES_IN_BUFFER;
postSP_data = bb.validSpikeBand(1 : SAMPLES_IN_BUFFER, :)';
clear bb

%%
postSP_normed = normalize(postSP_data, 'centered');

%% pwelch on normalized postSP data
Fs = 40000;

Pxx_postSP = [];
for nn = 1:32
    [ Pxx_postSP(:,nn), w ] = pwelch( postSP_normed( nn, : ), 40000*4 );
end

%%
%% plot
figure()
%ax(1) = subplot(2,1,1)
plot( w * Fs/(2*pi), mean(10 * log10( Pxx_postSP' )) );
xlabel( 'Freq ( Hz )' );
xlim([300 5000 ] );
title('Post_signalProcessing')










%% check peaks
figure
%plot( w * Fs/(2*pi), 10 * log10( Pxx_ICA(:, 5)' ) );
%xlim([0 5000])
%set(gcf, 'Position', [205 107 1377 859])
ch = 28:32;
for i = 1:5
    ax(i) = subplot(5,1,i);
    plot( w * Fs/(2*pi), 10 * log10( Pxx_ICA(:, ch(i))' ) );
    xlim([0 5000]);
    if i ~= 5
        set(gca, 'xtick', []);
    end
    title(['Component ' int2str(ch(i))])
end

xlabel( 'Freq ( Hz )' );
linkaxes(ax, 'xy')

%% peak detection
figure
x = w * Fs/(2*pi);
findpeaks(10 * log10(Pxx_normd(:, 1)'), x,  'MinPeakProminence', 5 )
xlim([0 5000])

%% verify on 5 random channels
figure
x = w * Fs/(2*pi);
ch = 1:5:25;
for i = 1:5
    ax(i) = subplot(5,1,i);
    findpeaks(10 * log10(Pxx_normd(:, ch(i))'), x,  'MinPeakProminence', 5 )
    xlim([0 5000]);
end
xlabel( 'Freq ( Hz )' );

linkaxes(ax, 'xy')

%% find peaks and apply notch filter
x = w * Fs/(2*pi);
for nn = 1:32
    [pks_tmp, locs_tmp] = findpeaks(10 * log10(Pxx_normd(:, nn)'), x,  'MinPeakProminence', 5 );
    pksInBand = find(locs_tmp <= 5000 & locs_tmp >= 150);
    locs{nn} = locs_tmp(pksInBand);
    pks{nn} = pks_tmp(pksInBand);
end

%%
dataFiltered = dataBuffer_normd;
parfor nn = 1:32
    for pk = 1:numel(locs{nn})
        freqToNotch = locs{nn}(pk); 
        W0 = freqToNotch/(Fs/2); BW = 2/(Fs/2); %BW = W0/35;
        [b,a] = iirnotch(W0, BW);
        dataFiltered(nn, :) = filter(b,a,dataFiltered(nn,:));
        disp(['finished ' int2str(pk)])
    end
end

%% find peaks on summed power spectrum across all channels and apply notch filter to each channel
x = w * Fs/(2*pi);
[pks, locs] = findpeaks(10 * log10(sum(Pxx_normd')), x,  'MinPeakProminence', 5 );
xlim([0 5000]);
pksInBand = find(locs <= 5000 & locs >= 300);
locs = locs(pksInBand);
pks = pks(pksInBand);

%% find peaks based on the summed data and apply notch filter
dataFiltered = dataBuffer_normd;
for pk = 1:numel(locs)
    freqToNotch = locs(pk); 
    W0 = freqToNotch/(Fs/2); BW = W0/35;
    [b,a] = iirnotch(W0, BW);
    for nn = 1:32
        dataFiltered(nn, :) = filter(b,a,dataFiltered(nn,:));
    end
    disp(['finished ' int2str(pk)])
end

%% pwelch on filtered data to verify
Fs = 40000;

Pxx_2 = [];
for nn = 1:32
    [ Pxx_2(:,nn), w ] = pwelch( dataBPFilted( nn, : ), 40000*4 );
end

%% compare to the unfiltered data
figure()
ax(1) = subplot(2,1,1)
plot( w * Fs/(2*pi), 10 * log10( sum(Pxx_normd') ) );
xlabel( 'Freq ( Hz )' );
xlim([300 5000 ] );
title('Before applying notch filter')

ax(2) = subplot(2,1,2)
plot( w * Fs/(2*pi), 10 * log10( sum(Pxx_2') ) );
xlabel( 'Freq ( Hz )' );
xlim([300 5000 ] );
title('After applying notch filter')

set(gcf, 'Position', [209 86 1288 863])
linkaxes(ax, 'xy')

%% compare to the unfiltered data
figure()
ax(1) = subplot(2,1,1)
plot( w * Fs/(2*pi), 10 * log10( Pxx_normd(:,2)' ) );
xlabel( 'Freq ( Hz )' );
xlim([0 7000 ] );
title('Before applying notch filter')

ax(2) = subplot(2,1,2)
plot( w * Fs/(2*pi), 10 * log10( Pxx_2(:, 2)' ) );
xlabel( 'Freq ( Hz )' );
xlim([0 7000 ] );
title('After applying notch filter')

set(gcf, 'Position', [209 86 1288 863])
linkaxes(ax, 'xy')



%% do bandpass filter and plot
filtLowCutoff = 300;
filtHighCutoff = 5000;
[b,a] = butter(4, [filtLowCutoff filtHighCutoff] / (Fs / 2), 'bandpass');
dataBPFilted = filtfilt(b, a, dataFiltered');
dataBPFilted = dataBPFilted';


%% verify before Pwelch
figure();
ax1 = subplot(2,1,1);
plot( dataBuffer(10,:) );
xlabel( 'time' );
ax2 = subplot(2,1,2);
plot( dataBuffer_normd(10,:) );
xlabel( 'time' );
linkaxes([ax1, ax2], 'xy');



%% verify summed Pxx
figure();
ax1 = subplot(2,1,1);
plot( w * Fs/(2*pi), 10 * log10( sum(Pxx') ) );
xlabel( 'Freq ( Hz )' );
xlim([0 2000 ] );
ax2 = subplot(2,1,2);
plot( w * Fs/(2*pi), 10 * log10( sum(Pxx_normd') ) );
xlabel( 'Freq ( Hz )' );
xlim([0 2000 ] );
linkaxes([ax1, ax2], 'xy');


%% verify single channel
figure();
ax1 = subplot(2,1,1);
plot( w * Fs/(2*pi), 10 * log10( Pxx(:, 1)' ) );
xlabel( 'Freq ( Hz )' );
xlim([0 5000 ] );
ax2 = subplot(2,1,2);
plot( w * Fs/(2*pi), 10 * log10( Pxx_normd(:, 1)' ) );
xlabel( 'Freq ( Hz )' );
xlim([0 5000 ] );
linkaxes([ax1, ax2], 'xy');



%%

xlim([0 5000 ] );
%% plot allchannel sum
figure()
plot( w * Fs/(2*pi), 10 * log10(  Pxx_normd(:,1)' ) );
xlabel( 'Freq ( Hz )' );
xlim([0 5000 ] );


%% plot individual channels
figure()
plot( w * Fs/2, 10 * log10( Pxx(:,1)') );
xlabel( 'Freq ( Hz )' );
xlim([0 2000 ] );

%% apply Reza's PLR code

dataFiltered = removePLI_multichan(dataBuffer, Fs, 27, [5,0.0001,7], [0.05,6,6], 0.7, 60);

%%
Pxx_2 = [];
for nn = 1:32
    [ Pxx_2(:,nn), w_2 ] = pwelch( dataFiltered( nn, fs*5:end ) , 40000*4);
end


%% notch filter
F0 = [188.3920, 302.0025, 376.7840, 503.3374, 565.6554, 604.0049, 611.1955, 654.3387];
W0 = (190/(95/30))/(Fs/2); BW = W0/35;
[b,a] = iirnotch(W0, BW);
for nn = 1:32
    dataFiltered(nn, :) = filter(b,a,dataBuffer(nn,:));
end

%%
Fs = 40000;

Pxx_2 = [];
for nn = 1:32
    [ Pxx_2(:,nn), w_2 ] = pwelch( dataFiltered( nn, : ), 40000*4 );
end

%%
figure()
subplot(2,1,1)
plot( w * Fs/(2*pi), 10 * log10( sum(Pxx') ) );
xlabel( 'Freq ( Hz )' );
xlim([0 5000 ] );
title('Before applying notch filter')

subplot(2,1,2)
plot( w_2 * Fs/(2*pi), 10 * log10( sum(Pxx_2') ) );
xlabel( 'Freq ( Hz )' );
xlim([0 5000 ] );
title('After applying notch filter')

set(gcf, 'Position', [209 86 1288 863])