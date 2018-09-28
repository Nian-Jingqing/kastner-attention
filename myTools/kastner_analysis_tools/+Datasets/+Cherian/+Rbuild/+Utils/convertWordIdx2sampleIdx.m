function [sample_idx] = convertWordIdx2sampleIdx(word_idx,words,sampleRate)
    time = words(word_idx,1);
    sample_idx = Datasets.Cherian.Rbuild.Utils.convert2sampleNumber(time,sampleRate);
end
