function spikeband_to_spikes_func(bb, date, multiple)
UE_baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/';
UE_file = ['UE_' date '.mat'];
load_UE_path = fullfile(UE_baseDir, UE_file);
clear UE
%load_UE_path = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/UE_03062019.mat';
load(load_UE_path);

startInds = UE.fixOn;
stopInds = UE.trialEnd + 400;
trialstruct = struct;
for itrial = 1:numel(startInds)
    trialstruct(itrial).isErrorTrial = UE.isErrorTrial(itrial);
    trialstruct(itrial).isEarlyError = UE.isEarlyError(itrial);
    trialstruct(itrial).cueType = UE.cueType(itrial);
    trialstruct(itrial).barType = UE.cueType(itrial);
    trialstruct(itrial).fixType = UE.fixType(itrial);
    trialstruct(itrial).isValidTarget = UE.isValidTarget(itrial);
    trialstruct(itrial).isSameObjTarget = UE.isSameObjTarget(itrial);
    trialstruct(itrial).isDiffObjTarget = UE.isDiffObjTarget(itrial);    
    trialstruct(itrial).startTime = startInds(itrial);
    trialstruct(itrial).endTime = stopInds(itrial);
    trialstruct(itrial).startInd = startInds(itrial);
    trialstruct(itrial).endInd = stopInds(itrial);
    trialstruct(itrial).condition = UE.cueType(itrial); % need to change!!        
end


%% get threshold to use for each channel
thresh = zeros(1, size(bb.minSpikeBand, 2)); % for each channel, there is one threshold
seg_len = 9062008;
num_segs = 5;
for iseg = 1:num_segs
    startInd = ((iseg - 1) * seg_len + 1);
    endInd = startInd + seg_len - 1;
    tmp = bb.validSpikeBand(startInd : endInd, :);
    wfstds = std(tmp, 0, 1);
    thresh = thresh + wfstds;
    clear tmp
    iseg
end

thresh = -multiple * (thresh / num_segs);

%% get var and remove NaN, also remove NaN for minSpikeBand
for ich = 1:size(bb.minSpikeBand,2)
    chVar{ ich } = bb.meanSquared( bb.meanSquaredChannel == ich );
    whereNan = find( isnan( chVar{ ich } ) );
    chVar{ ich } = chVar{ ich }( 1 : (whereNan(1) - 1) );
    %    whereNan_msb = find( isnan( bb.minSpikeBand_rmCA( : , ich ) ) );
    whereNan_msb = find( isnan( bb.minSpikeBand( : , ich ) ) );
    %chMsb{ ich } = bb.minSpikeBand_rmCA(1:(whereNan_msb(1) - 1), ich);
    chMsb{ ich } = bb.minSpikeBand(1:(whereNan_msb(1) - 1), ich);
end


%% get spikes
spikes = sparse(size(chMsb{1}, 1), numel(chMsb));
for ich = 1:size(bb.minSpikeBand,2)
    chThres = thresh(ich);
    leftValue = chMsb{ich} - chThres;
    spikes(leftValue <= 0, ich) = 1;
end
% remove spikes if show on more than 25 TCs (~80%)
allSpikesMS = sum(spikes,2);
spikes((allSpikesMS > 25), :) = 0;
% % get rid of MU 1 - 4
%spikes = spikes(:, 5:end);
stream.spikes = sparse(spikes);

%%
dtMS = 1;
C = Continuous.Continuous(stream, dtMS);
sigma_neural = 10;
C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );
r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );
%keyboard
%% save
savedir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/';

%cd(savedir);
saveFileName = ['LIP_spiking_r_', date, '.mat'];
outputfile = fullfile(savedir, saveFileName);
save(outputfile,'r', '-v7.3');
clear all