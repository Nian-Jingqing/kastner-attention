addpath('../Matlab_Offline_Files_SDK/');
addpath('rawDataProcessing');
%%
fn = 'Remy_RP_02182019_LIP_WB.pl2';
raw_data_path = '/mnt/scratch/feng/LIP';

%%

full_fn = fullfile( raw_data_path, fn );
pl2 = PL2GetFileIndex( full_fn );
%PL2Print(pl2.AnalogChannels);

%%
inputfile = full_fn;

% display information about the plx file
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information( inputfile );


disp(['Opened File Name: ' OpenedFileName]);
disp(['Version: ' num2str(Version)]);
disp(['Frequency : ' num2str(Freq)]);
disp(['Comment : ' Comment]);
disp(['Date/Time : ' DateTime]);
disp(['Duration : ' num2str(Duration)]);
disp(['Num Pts Per Wave : ' num2str(NPW)]);
disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);

% some information is only filled if the plx file version is >102
if ( Version > 102 )
    if ( Trodalness < 2 )
        disp('Data type : Single Electrode');
    elseif ( Trodalness == 2 )
        disp('Data type : Stereotrode');
    elseif ( Trodalness == 4 )
        disp('Data type : Tetrode');
    else
        disp('Data type : Unknown');
    end

    disp(['Spike Peak Voltage (mV) : ' num2str(SpikePeakV)]);
    disp(['Spike A/D Resolution (bits) : ' num2str(SpikeADResBits)]);
    disp(['Slow A/D Peak Voltage (mV) : ' num2str(SlowPeakV)]);
    disp(['Slow A/D Resolution (bits) : ' num2str(SlowADResBits)]);
end

% there's some other info available. I don't really know what it means.
tmp = PL2GetFileIndex( inputfile );



%%% GLOBALS
% samplerate is kept in the Freq variable
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

% half of the "prior" buffer is to deal with edge effects
% the other half is to provide valid data for the previous frame's edge effects
EDGE_EFFECT_BUFFER = SAMPLES_PER_MS * 1000;
PRIOR_DATA_BUFFER_SAMPLES = EDGE_EFFECT_BUFFER * 2;
POST_DATA_BUFFER_SAMPLES = EDGE_EFFECT_BUFFER;

%% get some more info about the file
% 'slowcounts' tells us how many analog samples exist for each channel
[tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);
% find channels that have analog samples defined
analogChannels = find(slowcounts);
numChannels = numel( analogChannels );

% how many analog samples are there?
numSamples = unique( slowcounts( analogChannels ) );

% this better be the same for all channels, otherwise we don't know how to process.
assert( numel( numSamples ) == 1, 'different samples have different lengths. don''t know how to process' );

disp( sprintf('Num Analog Samples: %g', numSamples ) );

% keep some amount of previous data around, for edge effects
% create a buffer for prior data
priorData = zeros( PRIOR_DATA_BUFFER_SAMPLES, ...
                   numChannels, 'single');

% create a buffer to store moving avg of mean square value
MOV_AVG_BUFFER_SIZE = 100;
prevMsVals = zeros(MOV_AVG_BUFFER_SIZE, numChannels, 'single');




%%
% initialize some variables for our big for loop
currentBufferStartInd = 1;
numFramesProcessed = 0;
channelMsInd = 1;
lastBlockTime = 0;

[ dataBuffer, samplesRead ] = readPL2Samples( inputfile, SAMPLES_IN_BUFFER, analogChannels, numFramesProcessed == 0 );

%% pwelch on raw dataBuffer
Fs = 40000;

Pxx = [];
for nn = 1:32
    [ Pxx(:,nn), w ] = pwelch( dataBuffer( nn, : ), 40000*4 );
end

%% normalize power
dataBuffer_normd = normalize(dataBuffer, 'centered');

%% perform ICA
dataBuffer_ICA = fastICA(dataBuffer_normd,32);

%% pwelch on normalized dataBuffer
Fs = 40000;

Pxx_normd = [];
for nn = 1:32
    [ Pxx_normd(:,nn), w ] = pwelch( dataBuffer_normd( nn, : ), 40000*4 );
end

%% pwelch on ICA dataBuffer
Fs = 40000;

Pxx_ICA = [];
for nn = 1:32
    [ Pxx_ICA(:,nn), w ] = pwelch( dataBuffer_ICA( nn, : ), 40000*4 );
end


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

