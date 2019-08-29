%filename = '/mnt/scratch/feng/Remy_02262019_PUL_Raw.pl2';
%outfile = '/mnt/scratch/cpandar/Remy_02262019_PUL_spikeband.mat';
%
%tic;
%broadband2streamMinMax( filename, outfile )
%toc;
%% first load UEs
% load UE and get trialstruct
clear UE 
load_UE_path = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/UE_05022019.mat';
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

%%
spikeBandFile = '/snel/share/share/data/kastner/Manoj/LIP/spikeBand/Remy_05022019_LIP_spikeband.mat';
bb = load(spikeBandFile);
bb = bb.spikeband;

%% get var and remove NaN, also remove NaN for minSpikeBand
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

%% remove common average
bb.minSpikeBand_rmCA = bb.minSpikeBand - mean(bb.minSpikeBand,2);

sd_factor = 1.5;
spikes = sparse(size(chMsb{1}, 1), numel(chMsb));
for ich = 1:size(bb.minSpikeBand,2)
    chMean = mean(chMsb{ich});
    chThres = chMean - sd_factor*chStd_cons{ich};
    leftValue = chMsb{ich} - chThres;
    spikes(leftValue <= 0, ich) = 1;
end
% remove spikes if show on more than 25 TCs (~80%)
allSpikesMS = sum(spikes,2);
spikes((allSpikesMS > 25), :) = 0;
% % get rid of MU 1 - 4
%spikes = spikes(:, 5:end);
stream.spikes = sparse(spikes);
dtMS = 1;
C = Continuous.Continuous(stream, dtMS);
sigma_neural = 10;
C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );
r = Datasets.PulvinarTools.pulvinarData( C.makeTrialsFromData( startInds, stopInds, trialstruct ) );
tc_r = R.Rstruct(r.r);

%%
saveDriftdir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/driftingCheck/05022019/';
if ~isdir(saveDriftdir)
    mkdir(saveDriftdir)
end
cd(saveDriftdir)

%% calculate mean FR in each trial
mean_FR = zeros(size(tc_r.r(1).spikes, 1), numel(tc_r.r));
for i = 1:numel(tc_r.r)
    mean_FR(:, i) = 1000*mean(tc_r.r(i).spikes, 2);
end


%% make plots
figure
x = 1:numel(tc_r.r);
for n = 1:size(tc_r.r(1).spikes, 1)
    %for n = 1
    clf;
    scatter(x, mean_FR(n,:))
    xlabel('Trials')
    ylabel('Mean FR')
    title(['Neuron ' int2str(n)])
    print(gcf,['Neuron ' int2str(n)], '-dpng');
end

