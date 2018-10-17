function [AvgFiringRate_InRF,AvgFiringRate_OffRF, SmoothedSpikingMatrix, TrialIndexInRF, TrialIndexOffRF, nTimesLFADS] = computeFRForPSTH(assembled, infield, sigma, binsize, olapChopped, type, smoothOrNot)
%COMPUTEFRFORPSTH Summary of this function goes here
%   Detailed explanation goes here

if strcmp(type, 'lfads')
    cutOff_LFADSRates = ceil(sigma/binsize);
    assembled.postSmoothCutOff( infield, 'rates_cutOff', cutOff_LFADSRates);
    SmoothedSpikingMatrix = permute(cat(3, assembled.r.rates_cutOff), [3 2 1]);
    nTimesLFADS = size(assembled.r(1).rates_cutOff, 2);
elseif strcmp(smoothOrNot, 'raw')
    cutOff_smoothedSpiking = binsize*ceil(sigma/binsize);
    assembled.postSmoothCutOff(infield, 'rawSpike_cutOff', cutOff_smoothedSpiking);
    % assembled.smoothFieldInR( infield, 'spike_smoothed', sigma, 1);
    %     assembled.postSmoothCutOff( 'spike_smoothed', 'spike_cutOff', cutOff_smoothedSpiking);
    rbinned = assembled.binData({'rawSpike_cutOff'}, [binsize]);
    SmoothedSpikingMatrix = permute(cat(3, rbinned.rawSpike_cutOff), [3 2 1]);
    nTimesLFADS = size(rbinned(1).rawSpike_cutOff,2);
else
    cutOff_smoothedSpiking = binsize*ceil(sigma/binsize);
    assembled.smoothFieldInR( infield, 'spike_smoothed', sigma, 1);
    assembled.postSmoothCutOff( 'spike_smoothed', 'smoothedSpike_cutOff', cutOff_smoothedSpiking);
    rbinned = assembled.binData({'smoothedSpike_cutOff'}, [binsize]);
    SmoothedSpikingMatrix = permute(cat(3, rbinned.smoothedSpike_cutOff), [3 2 1]);
    nTimesLFADS = size(rbinned(1).smoothedSpike_cutOff,2);
end

clVector = arrayfun(@(x) x.cueLoc, olapChopped);
rfLoc = olapChopped(1).rfloc;

AllCueLoc = unique(clVector);
rebinSize = binsize; % data was rebinned from 1ms to 10ms during LFADS

nTrials = length(assembled.r); % get trial number
nTimesRaw = size(assembled.r(1).(infield), 2); % get trial length for raw data, AKA, before re-binned
nNeurons = size(assembled.r(1).(infield), 1); % get neuron nubmer

AvgFiringRate = zeros(nNeurons, nTimesLFADS);
% initialize a avg firing rate matrix to store avg true firing rate (nNeurons x n Rebinned Times)
% for n = 1:nNeuron % loop over all neurons
% AvgLFADSRate = zeros(nNeurons, nTimesLFADS);
% initialize a avg LFADS rate matrix to store avg LFADS rates (nNeurons x n Rebinned Times)
AvgFiringRate_InRF = zeros(nNeurons, nTimesLFADS);
AvgFiringRate_OffRF = zeros(nNeurons, nTimesLFADS);

% % initialize avg Firing rates matrix for differet cue location
% AvgLFADSRate_InRF = zeros(nNeurons, nTimesLFADS);
% AvgLFADSRate_OffRF = zeros(nNeurons, nTimesLFADS);

% initialize avg LFADS rates matrix for differet cue location
TrialIndexInRF = false(nNeurons, nTrials);
TrialIndexOffRF = false(nNeurons, nTrials);
% initialize a matrix for each cue location to store the trial index. Ones
% would be index for trials that have the corresponding cue location

for n = 1:nNeurons % loop over all neurons
    if strcmp(type, 'lfads')
        AvgFiringRate(n,:) = (sum(SmoothedSpikingMatrix(:,:,n),1))*(1/nTrials);
    else
        AvgFiringRate(n,:) = (sum(SmoothedSpikingMatrix(:,:,n),1))*(1/nTrials)*(1000/rebinSize);
    end
    % Store the avg firing rate to the avgFiringRate matrix
%     AvgLFADSRate(n,:) = (sum(WholeLFADSRatesMatrix(:,:,n),1))*(1/nTrials);

    TrialIndexInRF(n,:) = clVector == rfLoc(n,1);
    nInRF = length(clVector(clVector == rfLoc(n,1)));
    % find the number of trials that cue loc is in RF
    TrialIndexOffRF(n,:) = clVector == rfLoc(n,2);
    % Find the index of the trials that cue loc is in Rf or off Rf
    nOffRF = length(clVector(clVector == rfLoc(n,2)));
    % find the number of trials that cue loc is Off RF
%     OutRF = AllCueLoc(~ismember(AllCueLoc, rfLoc(n,:))); 
%     % Find the cue locations that are not in or off Rf of this neuron
%     TrialIndexOutRF1(n,:) = clVector == OutRF(1);
%     % always find the index of the trials that cue loc number is smaller and
%     % put these index into the OutRF1
%     nOutRF1 = length(clVector(clVector == OutRF(1)));
%     % find the number of trials that cue loc is OutRF 1
%     TrialIndexOutRF2(n,:) = clVector == OutRF(2);
%     % always find the index of the trials that cue loc number is larger and
%     % put these index into the OutRF2
%     nOutRF2 = length(clVector(clVector == OutRF(2)));
%     % find the number of trials that cue loc is OutRF 2
    if strcmp(type, 'lfads')
        AvgFiringRate_InRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF);
        AvgFiringRate_OffRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF);
    else
        AvgFiringRate_InRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF)*(1000/rebinSize);
        AvgFiringRate_OffRF(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF)*(1000/rebinSize);
    end
%     AvgFiringRate_OutRF1(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1)*(1000/par.spikeBinMs);
%     AvgFiringRate_OutRF2(n,:) = (sum(SmoothedSpikingMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2)*(1000/par.spikeBinMs);

%     AvgLFADSRate_InRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexInRF(n,:),:,n),1))*(1/nInRF);
%     AvgLFADSRate_OffRF(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOffRF(n,:),:,n),1))*(1/nOffRF);
%     AvgLFADSRate_OutRF1(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF1(n,:),:,n),1))*(1/nOutRF1);
%     AvgLFADSRate_OutRF2(n,:) = (sum(WholeLFADSRatesMatrix(TrialIndexOutRF2(n,:),:,n),1))*(1/nOutRF2);

end

end

