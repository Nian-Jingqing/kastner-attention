function [cond_avg, channel_dayID] = getCondAvgAllDays(dataset, UE, alf, rateOrSpike, timePoints, nTimes, nConds, conditions, binsize_rescaled)
% cond_avg is the condition avg matrix (nChannel x nConditions*nTimes).
% There are 2 bar types and 8 cue types, so combined 16 conditions at cueOnset
% 16 is hard coded -- need to change
num_channels = zeros(1, numel(dataset));
for nday = 1:numel(dataset)
    num_channels(nday) = size(alf{nday}(1).(rateOrSpike), 1);
end
total_num_channels = sum(num_channels);

%cond_avg = zeros(length([good_channels{:}]), nConds*nTimes);
cond_avg = zeros(total_num_channels, nConds*nTimes);
% channel_dayID stores the day ID for each channel in the cond_avg matrix
%channel_dayID = zeros(1, length([good_channels{:}]));
channel_dayID = zeros(1, total_num_channels);

past_channelNum = 0; % how many channels in total in the past days
this_cond_where = 1:nTimes; % this is the inds of where the current cond is supposed to stored in the cond_avg matrix

% loop over each day and each bar type and each cue type
for nday = 1:numel(dataset)
    for nBar = 1:numel(conditions.barOn)        
        for nCue = 1 : numel(conditions.cueOn)
            % first select the condition to calculate 
            cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
            trialIndicesPerCond.cueOn.(cue_cond_field) = (UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue);

            % compute cond_avg of that condition. Used previous code when calculating PSTH
            [psth, raster_tensor, stderr] = preparePSTH_Manoj(alf{nday}, rateOrSpike, trialIndicesPerCond.cueOn.(cue_cond_field), [alf{nday}.cueOnset], timePoints, binsize_rescaled);
            cond_avg((1:num_channels(nday)) + past_channelNum, this_cond_where) = psth;
            this_cond_where = this_cond_where + nTimes;
            % update this_cond_where to indicate the inds for the next condition
        end
    end
    % fill in channel_dayID for the channels of the current day
    channel_dayID((1:num_channels(nday)) + past_channelNum) = nday;
    past_channelNum = past_channelNum + num_channels(nday);
    this_cond_where = 1:nTimes;
end