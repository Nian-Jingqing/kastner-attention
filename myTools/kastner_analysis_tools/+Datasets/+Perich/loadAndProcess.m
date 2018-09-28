function [C , startInds, stopInds, trialstruct ]= Perich_loadAndProcess ( fname )

% load matt's data
d = load( fname);
d = d.data;


%% convert kinematic data into continuous data stream

% time in bins
total_samples = size(d.kin.t, 1);
% time in seconds
total_time = d.meta.dataWindow(2);

%samples_per_s = total_samples / total_time;
samples_per_s = round( total_samples / total_time );
desired_samplerate = 1000; % (Hz)

% i don't yet understand the d.kin.t field
apparent_samplerate = mean(diff(d.kin.t));


% iterate over all kinematic fields and resample them
res_kin = [];
% extract the fields in in the kinematic data struct 
fs = fields(d.kin);
for nf = 1:numel( fs )
    res_kin.( fs{nf} ) = resample( d.kin.( fs{nf} ), desired_samplerate, ...
                                         samples_per_s );
end

%
% destination for the big stream of data:
stream = [];

% turn kinematics into one big array
%    format: [ x y vx vy ax ay ]
stream.kin = [res_kin.x(:) ...
              res_kin.y(:) ...
              res_kin.vx(:) ...
              res_kin.vy(:) ...
              res_kin.ax(:) ...
              res_kin.ay(:)];

% turn spiketimes into a (sparse) array of 1s and 0s
stream.spikes = sparse(total_samples, numel(d.units) );
for iunit = 1:numel(d.units)
    % get the spike times
    spks = d.units( iunit ).spikes.ts;
    % trim to times within the data window
    spksshort = spks(spks >= d.meta.dataWindow(1) & ...
                     spks < d.meta.dataWindow(2) );
    spk_time_inds = 1 + floor( spksshort * desired_samplerate );
    stream.spikes(spk_time_inds, iunit) = 1;
end


%% turn trial table into a struct array
numTrials = size( d.trials, 1 );
trialstruct = struct;
startInds = zeros(numTrials, 1);
stopInds = zeros(numTrials, 1);
for itrial = 1 : numTrials
    % these times are in seconds
    trialStart = d.trials.startTime(itrial);
    trialStop = d.trials.endTime(itrial);

    startInds(itrial) = round( trialStart * desired_samplerate );
    stopInds(itrial) = round( trialStop * desired_samplerate );

    % store these times in the trial struct
    trialstruct(itrial).startTime = trialStart;
    trialstruct(itrial).endTime = trialStop;
    trialstruct(itrial).startInd = startInds(itrial);
    trialstruct(itrial).endInd = stopInds(itrial);
    
    % trial number
    trialstruct(itrial).trialNum = d.trials.number(itrial);
    % trial result
    trialstruct(itrial).result = d.trials.result(itrial);
    % conditionID
    trialstruct(itrial).conditionID = d.trials.tgtID(itrial);
    % target position
    trialstruct(itrial).posTarget = d.trials.tgtCtr(itrial)';
    % target angle
    trialstruct(itrial).angleTarget =  d.trials.tgtDir(itrial);

    % get event times relative to trialStart (in seconds), then convert to samples
    % target onset
    targetOnset = d.trials.tgtOnTime(itrial)- trialStart;
    trialstruct(itrial).targetOnset = round( targetOnset * desired_samplerate );
    % go cue
    goCue = d.trials.goCueTime(itrial) - trialStart;
    trialstruct(itrial).goCue = round( goCue * desired_samplerate );

    % target corners is also a thing, but I don't really care about it.
end
    
    
    


% use "continuous" class to manipulate data from here on out
dtMS = 1; % data is sampled at 1 ms
C = Continuous.Continuous( stream , dtMS);

