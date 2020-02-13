%%


%dataset(1).date = '02272019';
%dataset(2).date = '03042019';
%dataset(3).date = '03072019';
%dataset(4).date = '03092019';
%dataset(5).date = '03102019';
%dataset(6).date = '03122019';
%dataset(7).date = '03132019';
%dataset(8).date = '03152019';
dataset(1).date = '03082019';
dataset(2).date = '03272019';



addpath('../Matlab_Offline_Files_SDK/');
addpath('rawDataProcessing');
%%
day = 2;
fn = ['Remy_RP_' dataset(day).date '_LIP_WB.pl2'];
%fn = 'Remy_RP_03092019_LIP_WB.pl2';
%raw_data_path = '/mnt/scratch/feng/LIP';
raw_data_path = '/home/feng';
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/moreNewSessions/spectrogram/';
%raw_data_path = '/home/feng/';

full_fn = fullfile( raw_data_path, fn );
pl2 = PL2GetFileIndex( full_fn );
%PL2Print(pl2.AnalogChannels);
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
numFramesProcessed = 0; % was 0.
channelMsInd = 1;
lastBlockTime = 0;

while numFramesProcessed < 4
    [ dataBuffer, samplesRead ] = readPL2Samples( inputfile, SAMPLES_IN_BUFFER, analogChannels-128, numFramesProcessed == 0 );
    
    % first normalize
    dataBuffer_normd = normalize(dataBuffer, 'centered');

    % get Pxx for un-processed data
    Fs = 40000;
    Pxx_normd = [];
    for nn = 1:size(dataBuffer_normd, 1)
        [ Pxx_normd(:,nn), w ] = pwelch( dataBuffer_normd( nn, : ), 40000*4 );
    end
    x = w * Fs/(2*pi);

    % do lowpass - notch - bandpass
    %dataFiltered = lowPassNotchBandPass(dataBuffer_normd);

    % do notch - bandpass
    dataFiltered = notchThenBandPass(dataBuffer_normd);
    

    % get Pxx for processed data
    Pxx_post = [];
    for nn = 1:size(dataFiltered, 1)
        [ Pxx_post(:,nn), w ] = pwelch( dataFiltered( nn, : ), 40000*4 );
    end
    x = w * Fs/(2*pi);

    % generate comparison plot
    clf
    subplot(2,1,1)
    plot(x, 10*log10(Pxx_normd(:,1)'))
    %plot(x, mean(10*log10(Pxx_normd')))
    xlabel( 'Freq (Hz)');
    xlim([0 7000]);
    xticks([300]);
    xticklabels({'300'})
    ylim([-90 30])
    title('Pre SP')

    subplot(2,1,2)
    plot(x, 10*log10(Pxx_post(:,1)'))
    %plot(x, mean(10*log10(Pxx_normd_2')))
    xlabel( 'Freq (Hz)');
    xlim([0 7000]);
    xticks([300 1000 2000 3000 4000 5000 6000]);
    xticklabels({'300', '1000', '2000', '3000', '4000', '5000', '6000'})
    ylim([-90 30])
    title('Post SP')
    saveFileName = ['Day' dataset(day).date 'DataFrame' int2str(numFramesProcessed+1)]
    suptitle(saveFileName)

    cd(saveDir)
    print(gcf, saveFileName, '-dpng');
    
    numFramesProcessed = numFramesProcessed + 1;
    numFramesProcessed
end


