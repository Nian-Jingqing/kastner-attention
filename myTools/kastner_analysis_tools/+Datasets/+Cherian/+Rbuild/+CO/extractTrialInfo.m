function [ trialStruct ] = extractTrialInfo( words, spikeInfo, targStruct, stream, task, sampleRate )
% extracts trial information from madStruct and put in trial structure

% words : cue words associated with the experiment from madstruct
% targStruct : structure holding target info extracted from targ table
% stream : continous data that was extracted previous to this function
% sampleRate : sample rate of analog data

% trialStruct: struct holding relevant trial information for experiment
%addpath( '/snel/home/lwimala/Projects/MATLAB/CHERIAN/CODE-LW/')
    outputDataSize = size( stream.kin, 1 );
    % extract channels from words
    time = words( :, 1 );
    cue = words( :, 2 );

    % find start words
    cueTargOn = find( cue >= 64 & cue <= 71 );
    
    % If there are more target on cues than targets trim cue words
    if size( cueTargOn, 1 ) > size ( targStruct, 2 )
        cueTargOn = cueTargOn( 1 : size( targStruct, 2 ) );
    end
    
    % redundant assert that number of start words equal number of trials
    assert( size( cueTargOn, 1 ) == size( targStruct, 2 ),'ERROR: Number of target on words does not match number of trials.' )

    % intialize for loop
    trialStruct = struct();
    for i = 1:size( cueTargOn, 1 )
        
        % extract trial start times for current and next trial
        timeTrialStart = time( cueTargOn( i ) );
        assert( timeTrialStart == targStruct( i ).targOnTime, ...
                'ERROR: targOn  word time does not match targStruct time.')
        
        timeTrialGoCue = time( cueTargOn( i ) + 1 );
        assert( timeTrialGoCue == targStruct( i ).goCueTime, ...
                'ERROR: goCue word time does not match targStruct time.')
        
        if i == size(cueTargOn,1)
            idxNextTrialStart = outputDataSize + 1;
        else
            timeNextTrialStart = time( cueTargOn( i + 1 ) );
            idxNextTrialStart = Rbuild.Utils.convert2sampleNumber( timeNextTrialStart, sampleRate );
        end
        
        % convert times to sample indices
        idxTrialStart = Rbuild.Utils.convert2sampleNumber( timeTrialStart, sampleRate );
        idxTrialGoCue = Rbuild.Utils.convert2sampleNumber( timeTrialGoCue, sampleRate );
        idxTrialEnd = idxNextTrialStart - 1;

        % extract trial cue words for target onset and go
        cueTrialTargOn = cue( cueTargOn( i ) );
        cueTrialGo = cue( cueTargOn( i ) + 1 );

        % check target onset to determine condition ID
        switch cueTrialTargOn
          case 64
            trialCondId = 1;
          case 65
            trialCondId = 2;
          case 66
            trialCondId = 3;
          case 67
            trialCondId = 4;
          case 68
            trialCondId = 5;
          case 69
            trialCondId = 6;
          case 70
            trialCondId = 7;
          case 71
            trialCondId = 8;
        end

        % check go cue to determine if trial is success or failure
        switch cueTrialGo
          case 49
            trialType = TrialType.Success;
          case 33
            trialType = TrialType.Failure;
        end
        trialStruct( i ).trialNum = i;
        trialStruct( i ).trialType = trialType;
        trialStruct( i ).targOnTime = targStruct( i ).targOnTime;
        trialStruct( i ).goCueTime = targStruct( i ).goCueTime;
        trialStruct( i ).targPos = targStruct( i ).pos;
        trialStruct( i ).startIdx = idxTrialStart;
        trialStruct( i ).endIdx = idxTrialEnd;
        trialStruct( i ).absTargOn = idxTrialStart;
        trialStruct( i ).absGoCue = idxTrialGoCue;
        trialStruct( i ).relTargOn = 1;
        trialStruct( i ).relGoCue = idxTrialGoCue - idxTrialStart;
        trialStruct( i ).conditionId = trialCondId;
        trialStruct( i ).task = task.Task;
        trialStruct( i ).spikeInfo = spikeInfo;
        trialStruct( i ).sampleRate = sampleRate;
        try
            trialStruct( i ).force = task.force;
        catch ME
            trialStruct( i ).force = 'None';
        end
    end
end
