function [avg_rates, raster_tensor, stderror] = preparePSTH_Manoj(day_struct, infield, trial_indices, event_times, timePoints, rebinSize)
% This function compute averaged firing rate across trials for each neuron, as well as organize the data into a tensor (n x T x time)
%   Detailed explanation goes here
% this function takes in a structure where each element is a trial, with the
% fields of LFADS rates, LFADS factors, smoothed spikes, and key event times
% for each trial. This function also takes in which field we want to prepare
% for PSTH. This function also takes in trial indices we want for the PSTH
% and around what event time we want for the PSTH. it also takes in a
% "timePoints" which specify the time around the event. This function will return
% a matrix (n x t), each n is the averaged firing rate of the neuron across
% trials. This function also returns a tensor (n x T x time) for the raster.

%halfWide = round(howWide/2);
%window = round([-1*halfWide  halfWide] / binsize_rescaled);
%timePoints = window(1):window(2);
keepTrials_struct = day_struct(trial_indices);
event_times = event_times(trial_indices);
nNeurons = size(keepTrials_struct(1).(infield), 1);
avg_rates = zeros(nNeurons, length(timePoints));
raster_tensor = zeros(nNeurons, numel(keepTrials_struct), length(timePoints));
for n = 1:nNeurons
    for itrial = 1:numel(keepTrials_struct)
        raster_tensor(n, itrial, :) = keepTrials_struct(itrial).(infield)(n, event_times(itrial) + timePoints);
    end
end
if strcmp(infield, 'rates') || strcmp(infield, 'factors')
    avg_rates = squeeze(mean(raster_tensor, 2));
    tmp_std = squeeze(std(raster_tensor,0, 2));
    stderror = tmp_std/sqrt(numel(keepTrials_struct));
else
    avg_rates = squeeze(mean(raster_tensor, 2))*(1000/rebinSize);
    tmp_std = squeeze(std((1000/rebinSize)*raster_tensor,0, 2));
    stderror = tmp_std/sqrt(numel(keepTrials_struct));
    raster_tensor = raster_tensor*(1000/rebinSize);
end

end
