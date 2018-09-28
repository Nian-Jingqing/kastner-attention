function [words] = extractWords( mad, mA, sampleRate )
    
    outputDataSize = numel(mA(1).waveform);
    % get analog start time
    try 
        analogStartTime = mA(1).info.starttime;
    catch ME
        analogStartTime = mad.analoginfo(1).starttime;
    end
    % get words
    try 
        tmp_words = mA(1).info.comment.words;
    catch ME
        tmp_words = mad.analoginfo( 1 ).comment.words;
    end
    
        
    time = (1:outputDataSize)/sampleRate;

    tmp_words(:,1) = tmp_words(:,1) - analogStartTime;
    % Set "start" word at beginning of data so trial table will include all data
    words = tmp_words;
end
