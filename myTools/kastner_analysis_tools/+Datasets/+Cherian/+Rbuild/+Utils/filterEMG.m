function [ abs_emg ] = filterEMG(emg,freqcutoffHigh, freqcutoffLow, sampleRate)

    % 4-Pole 20 Hz High Pass Butterworth Filter
    
    [BH,AH] = butter(4,freqcutoffHigh / (sampleRate/2),'high'); % CREATE High Pass Filter
    hp_emg = filter(BH,AH,emg');                               % APPLY High Pass Filter to EMG
    hp_emg = hp_emg';                                          % TRANSPOSE EMG channels to be in rows
    
    abs_emg = abs(hp_emg);                                     % RECTIFY EMG signals
    
    %    % 4-Pole 10 Hz Low Pass Butterworth Filter
    
    %    [BL,AL] = butter(4,freqcutoffLow / (sampleRate/2),'low');    % CREATE Low Pass Filter
    %    filt_emg = filter(BL,AL,abs_emg');                         % APPLY Low Pass Filter EMG channels
    %    filt_emg = filt_emg';                                      % TRANSPOSE EMG channels to be in rows

end

