function [sampleNumber] = convert2sampleNumber(sampleTime,sampleRate)
    sampleNumber = round(sampleTime.*sampleRate);
end
