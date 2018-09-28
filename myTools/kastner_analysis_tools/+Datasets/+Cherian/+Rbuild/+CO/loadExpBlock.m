function [ C, trialStruct, startInds, stopInds, exp_name ] = loadExpBlock( data , blockIdx )

% load in experiment block
expBlock = data( 1 ).Experiment_blocks;
block = expBlock{ blockIdx };

% extract madStruct and target table
mad  = block.madstruct;
targ = block.Target_table;
task = block.Task_info;

% Extract analog data and sample rate
mA = mad.analog;
sampleRate = mA( 1 ).info.samplerate;
% set new sample rate to 1 ms time bin
newSR = 1000;

% get experiment name
exp_name = strrep( mad.header.matfilename, '.mat', '' );

% Set relevant dataset parameters
excludeN = [];

offset_Xpos = -10;
offset_Ypos = 33.5;

% Extract experiment words
words = Rbuild.Utils.extractWords( mad, mA, sampleRate );

% Extract kinematic data
x = Rbuild.Utils.extractKin( mA, offset_Xpos, offset_Ypos );
x = Rbuild.Utils.resampleChannel( x, sampleRate, newSR );
% Extract raw EMG
tmp_emg = Rbuild.Utils.extractEMG( mA );

% Set filter cutoffs
freqCutHigh = 20;
freqCutLow = 10;

% Filter EMG 
emg = Rbuild.Utils.filterEMG( tmp_emg, freqCutHigh, freqCutLow, sampleRate );
emg = Rbuild.Utils.resampleChannel( emg, sampleRate, newSR );
% Clear raw emg 
clear tmp_emg

% Extract spikes
%[ spikes, spikeInfo ] = Rbuild.Utils.extractSpikes( mad, excludeN, sampleRate );
%spikes = Rbuild.Utils.rebinChannel( spikes, sampleRate, newSR );
%data = Rbuild.Utils.extractSpikes_CP( mad );
%data = Rbuild.Utils.resampleSpikes( mad, data, newSR );
data = Rbuild.Utils.extract_and_resample_spikes( mad, newSR );
%spikes = Rbuild.Utils.rebinChannel( data.spikes, sampleRate, newSR );
spikes = data.spikes;
spikeInfo = data.spikeInfo;
if size( spikes, 1 ) < size( emg, 1 )
    dataLength = size( spikes, 1 );
    emg = emg( 1 : dataLength, : );
    x = x( 1 : dataLength, : );
elseif size( emg, 1 ) > size( spikes, 1 )
    dataLength = size( emg, 1 );
    spikes = spikes( 1 : dataLength, : );
end

sampleRate = newSR;

% 
stream = [];

stream.kin = x;
stream.spikes = spikes;
stream.emg = emg;

[ targStruct ] = Rbuild.CO.extractTargetInfo( targ, mA);

[ trialStruct ] = Rbuild.CO.extractTrialInfo( words, spikeInfo, targStruct, stream, ...
                                              task, sampleRate );


dtMS = (1/sampleRate)*1000;

startInds = [ trialStruct.startIdx ];
stopInds = [ trialStruct.endIdx ];

C = Continuous.Continuous(stream, dtMS);
    

