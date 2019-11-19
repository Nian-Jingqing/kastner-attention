%% set up day information
datasets(1).shortName = '0218';
datasets(1).midName = '021819';
datasets(1).longName = '02182019';
datasets(2).shortName = '0306';
datasets(2).midName = '030619';
datasets(2).longName = '03062019';
datasets(3).shortName = '0311';
datasets(3).midName = '031119';
datasets(3).longName = '03112019';
datasets(4).shortName = '0314';
datasets(4).midName = '031419';
datasets(4).longName = '03142019';
datasets(5).shortName = '0406';
datasets(5).midName = '040619';
datasets(5).longName = '04062019';
datasets(6).shortName = '0425';
datasets(6).midName = '042519';
datasets(6).longName = '04252019';
datasets(7).shortName = '0502';
datasets(7).midName = '050219';
datasets(7).longName = '05022019';

%% set up the factor of sd for threshold for each day
sd(1).factor = 1.75;
sd(2).factor = 1.75;
sd(3).factor = 1.75;
sd(4).factor = 1.75;
sd(5).factor = 1.75;
sd(6).factor = 1.5;
sd(7).factor = 1.75;

%% set up directories
spikeBandPath = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/';
load_UE_path = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/';
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/withoutDataMasking/higherThreshold/';

%%
for nday = 1:numel(datasets)
    spikeBandFileName = ['Remy_' datasets(nday).longName '_LIP_spikeband.mat'];
    spikeBandFile = fullfile(spikeBandPath, spikeBandFileName);
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
    % get spikes into stream
    clear stream
    stream.spikes = sparse(spikes);

    %% load pre-saved UE
    clear UE
    disp( sprintf( 'loading UEs %g / %g', 1, numel( datasets ) ) );
    fname = sprintf( '%sUE_%s.mat', load_UE_path, datasets( nday ).longName );
    load( fname );

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
    r = Datasets.PulvinarTools.pulvinarData_noExtInp( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );

    combinedData.r = r;
    combinedData.UE = UE
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    cd(saveDir);
    saveName = [datasets(nday).midName '_v1.mat'];
    %saveName = '020819_v1.mat';
    save(saveName, 'combinedData');
end





