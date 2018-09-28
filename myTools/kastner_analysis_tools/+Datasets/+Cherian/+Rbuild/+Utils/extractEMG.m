function [ emg ] = extractEMG( mA )
% extracts emg analog data from madStruct
    
% mA : mad.analog ( Analog field of madStruct )    
% emg : Returns EMG signals in column matrix
    
% flexor carpi radialis
    flCR = mA( 9 ).waveform;
    % medial biceps
    mBic = mA( 10 ).waveform;
    % lateral biceps
    lBic = mA( 11 ).waveform;
    % brachioradialis
    brac = mA( 12 ).waveform;
    % clavicular head of pectoralis
    cPec = mA( 13 ).waveform;
    % anterior deltoid
    aDel = mA( 14 ).waveform;
    % extensor carpi radialis
    exCR = mA( 15 ).waveform;
    % lateral triceps
    lTri = mA( 16 ).waveform;
    % latissimus dorsi
    latD = mA( 17 ).waveform;
    % medial deltoid
    mDel = mA( 18 ).waveform;

    emg = [ cPec latD aDel mDel mBic lBic lTri brac exCR flCR ];
end
