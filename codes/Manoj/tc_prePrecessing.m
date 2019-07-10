%% set up day information
datasets(1).shortName = '0208';
datasets(1).midName = '020819';
datasets(1).longName = '02082019';
datasets(2).shortName = '0218';
datasets(2).midName = '021819';
datasets(2).longName = '02182019';
datasets(3).shortName = '0226';
datasets(3).midName = '022619';
datasets(3).longName = '02262019';
datasets(4).shortName = '0227';
datasets(4).midName = '022719';
datasets(4).longName = '02272019';
datasets(5).shortName = '0308';
datasets(5).midName = '030819';
datasets(5).longName = '03082019';
datasets(6).shortName = '0310';
datasets(6).midName = '031019';
datasets(6).longName = '03102019';
datasets(7).shortName = '0311';
datasets(7).midName = '031119';
datasets(7).longName = '03112019';

%% set up the factor of sd for threshold for each day
sd(1).factor = 1.8;
sd(2).factor = 1.5;
sd(3).factor = 2;
sd(4).factor = 1.5;
sd(5).factor = 2;
sd(6).factor = 2;
sd(7).factor = 2;

%% set up directories
spikeBandDir = '/snel/share/share/data/kastner/Manoj/PUL/spikeBand/';
loadpath = ['/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/'];
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/tc/';

%%
for nday = 1:numel(datasets)
    spikeBandFileName = ['Remy_' datasets(nday).longName '_PUL_spikeband.mat'];
    spikeBandFile = fullfile(spikeBandDir, spikeBandFileName);
    bb = load(spikeBandFile);
    bb = bb.spikeband;
    for ich = 1:size(bb.minSpikeBand,2)
        chVar{ ich } = bb.meanSquared( bb.meanSquaredChannel == ich );
        whereNan = find( isnan( chVar{ ich } ) );
        chVar{ ich } = chVar{ ich }( 1 : (whereNan(1) - 1) );
        whereNan_msb = find( isnan( bb.minSpikeBand( : , ich ) ) );
        chMsb{ ich } = bb.minSpikeBand(1:(whereNan_msb(1) - 1), ich);
    end
    
    %% compute constant std
    chStdVec = zeros(1,32);
    for ich = 1:numel( chVar )
        %use mean of variance across entire session for each channel
        chStd_cons{ich} = sqrt( mean( chVar{ich} ) );
        chStdVec(ich) = chStd_cons{ich};
    end

    %% get spikes
    spikes = sparse(size(chMsb{1}, 1), numel(chMsb));
    for ich = 1:size(bb.minSpikeBand,2)
        chMean = mean(chMsb{ich});
        chThres = chMean - sd(nday).factor*chStd_cons{ich};
        leftValue = chMsb{ich} - chThres;
        spikes(leftValue <= 0, ich) = 1;
    end

    % remove spikes if show on more than 25 TCs (~80%)
    allSpikesMS = sum(spikes,2);
    spikes((allSpikesMS > 25), :) = 0;

    % get rid of MU 1 - 4
    spikes = spikes(:, 5:end);
    stream.spikes = sparse(spikes);

    %% load pre-saved data to get UE and external inputs
    disp( sprintf( 'loading chopped day %g / %g', 1, numel( datasets ) ) );
    fname = sprintf( '%s%s_v1.mat', loadpath, datasets( nday ).midName );
    tmp = load( fname );
    olapChopped_SS = tmp.combinedData;
    UE = tmp.combinedData.UE;

    %% trialize the threshold crossing (TC) data
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

    %% put spike train into a Continuous class
    dtMS = 1;
    C = Continuous.Continuous(stream, dtMS);
    sigma_neural = 10;
    C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );

    %% turn into a trialized (R) struct
    r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

    %% combine and save
    olapChopped_TC.r = r;
    for i = 1:numel(olapChopped_TC.r.r)
        olapChopped_TC.r.r(i).externalInputs = olapChopped_SS.r.r(i).externalInputs;
    end
    combinedData = olapChopped_TC;
    combinedData.UE = UE
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    cd(saveDir);
    saveName = [datasets(nday).shortName '19_v2.mat'];
    %saveName = '020819_v1.mat';
    save(saveName, 'combinedData');
end




