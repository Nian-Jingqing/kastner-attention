function [ C, trialStruct, startInds, stopInds, exp_name ] = loadExpBlock( data , blockIdx )

% load in experiment block
expBlock = data( 1 ).Experiment_blocks;
block = expBlock{ blockIdx };

% extract madStruct and target table
mad  = block.madstruct;
targ = block.Target_table;

% Extract analog data and sample rate
mA = mad.analog;
try
    sampleRate = mA( 1 ).info.samplerate;
catch ME
    sampleRate = mad.analoginfo(1).samplerate;
end

% set new sample rate to 1 ms time bin
newSR = 1000;
% Set relevant dataset parameters
excludeN = [];

offset_Xpos = -10;
offset_Ypos = 33.5;

% get experiment name
exp_name = strrep( mad.header.matfilename, '.mat', '' );

% Extract experiment words
words = Datasets.Cherian.Rbuild.Utils.extractWords( mad, mA, sampleRate );
words = [ [ 0.2500 double( 18 ) ]; words];
words = [ words; [ ( numel( mA( 1 ).waveform ) ) / sampleRate double( 244 ) ] ];
words = [ words; [ ( numel( mA( 1 ).waveform ) + 1 ) / sampleRate double( 18 ) ] ];

% Extract kinematic data
x = Datasets.Cherian.Rbuild.Utils.extractKin( mA, offset_Xpos, offset_Ypos );
x = Datasets.Cherian.Rbuild.Utils.resampleChannel( x, sampleRate, newSR );
% Extract raw EMG
tmp_emg = Datasets.Cherian.Rbuild.Utils.extractEMG( mA );

% Set filter cutoffs
freqCutHigh = 20;
freqCutLow = 10;

% Filter EMG 
emg = Datasets.Cherian.Rbuild.Utils.filterEMG( tmp_emg, freqCutHigh, freqCutLow, sampleRate );
emg = Datasets.Cherian.Rbuild.Utils.resampleChannel( emg, sampleRate, newSR );
% Clear raw emg 
clear tmp_emg

% Extract spikes
data = Datasets.Cherian.Rbuild.Utils.extract_and_resample_spikes( mad, newSR );
spikes = data.spikes;
spikeInfo = data.spikeInfo;

if size( spikes, 1 ) < size( emg, 1 )
    dataLength = size( spikes, 1 );
    emg = emg( 1 : dataLength, : );
    x = x( 1 : dataLength, : );
else
    dataLength = size( emg, 1 );
    spikes = spikes( 1 : dataLength, : );
end



if sum(spikeInfo( 23, : )  == [ 35 0 ]) == 2
    disp( 'Anil you bastard! Get those useless neurons OUT OF HERE!!!' )
    spikes( :, 23 ) = [];
    spikeInfo( 23, : ) = [];
end

% set new sample rate
sampleRate = newSR;

% 
stream = [];

stream.kin = x;
stream.spikes = spikes;
stream.emg = emg;

[ targStruct ] = Datasets.Cherian.Rbuild.FF.extract_target_info( targ, words, mA, mad );

sampleIndices = Datasets.Cherian.Rbuild.FF.generateTrialSampleIndices( targStruct, mA, words, sampleRate );

[ trialStruct ] = Datasets.Cherian.Rbuild.FF.extractTrialInfo( sampleIndices, spikeInfo, targStruct, sampleRate, mad );


dtMS = ( 1/sampleRate )*1000;

startInds = [ trialStruct.abstSidx ];
stopInds = [ trialStruct.abstEidx ];

C = Continuous.Continuous( stream, dtMS);
    

