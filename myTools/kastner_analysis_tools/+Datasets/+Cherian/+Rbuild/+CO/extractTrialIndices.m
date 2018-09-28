function [ trialStruct ] = extractTrialIndices( words, targStruct, stream, sampleRate )
    outputDataSize = size(stream.kin,2);
% Extract channels from words
    time = words( :, 1 );
    cue = words( :, 2 );

    % Find start words
    cueTargOn = find( cue >= 64 & cue <= 71 );

    % Assert that number of start words equal number of trials
    assert( size( cueTargOn, 1 ) == size( targStruct, 2 ),'ERROR: Number of target on words does not match number of trials.' )

    % Intialize for loop
    trialStruct = struct();
    for i = 1:size( cueTargOn, 1 )
        
        % Extract trial start times for current and next trial
        timeTrialStart = time( cueTargOn( i ) );
        timeTrialGoCue = time( cueTargOn( i ) + 1 );
        
        if i == size(cueTargOn,1)
            idxNextTrialStart = outputDataSize + 1;
        else
            timeNextTrialStart = time( cueTargOn( i + 1 ) );
            idxNextTrialStart = Rbuild.Utils.convert2sampleNumber( timeNextTrialStart, sampleRate );
        end
        
        % Convert times to sample indices
        idxTrialStart = Rbuild.Utils.convert2sampleNumber( timeTrialStart, sampleRate );
        idxTrialGoCue = Rbuild.Utils.convert2sampleNumber( timeTrialGoCue, sampleRate );
        idxTrialEnd = idxNextTrialStart - 1;

        % Extract trial cue words for target onset and go
        cueTrialTargOn = cue( cueTargOn( i ) );
        cueTrialGo = cue( cueTargOn( i ) + 1 );

        % Check target onset to determine condition ID
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

        % Check go cue to determine if trial is success or failure
        switch cueTrialGo
          case 49
            trialType = TrialType.Success;
          case 33
            trialType = TrialType.Failure;
        end

        trialStruct( i ).startIdx = idxTrialStart;
        trialStruct( i ).endIdx = idxTrialEnd;
        trialStruct( i ).trialType = trialType;
        trialStruct( i ).goCue = idxTrialGoCue;
        trialStruct( i ).condId = trialCondId;
    end
end
