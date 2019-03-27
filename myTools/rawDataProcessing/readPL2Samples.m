function [ dataBuffer, samplesRead ] = readPL2Samples( filename, samplesToRead, whichChannels, firstRead )
%% function [ dataBuffer, samplesRead ] = readPL2Samples( filename, samplesToRead, whichChannels, firstRead )
% Reads multichannel data from a PL2 file
% uses PL2 low-level functions for speed
%
% Function is called multiple times for the same file
% uses a persistent variable to pick up reading where it left off
%
% must pass in all arguments:
%   filename: PL2 filename
%   samplesToRead: How far into the file to read (relative to last call)
%   whichChannels: What channels do we expect to see (channel numbers)
%      If there are other channels in the dataset, you'll see errors
%   firstRead: Is this the first call to the function?
%      If you have called the same function for a different file,
%      this also provides a way to reset.
%
% 2019-03-27: initial version (CP)



% keep 'blockInfo' around for repeated calls to this function
persistent blockInfo;

% PL2 data block size
% obv hardcoded, but if this changes you'll notice some warning messages
NUM_SAMPLES_PLEXON_BLOCK = 65535;

% how many channels are there?
numChannels = numel( whichChannels );

dataBuffer = nan( numChannels, samplesToRead );

ich = 1;
samplesRead = 0;
expectedSamples = NUM_SAMPLES_PLEXON_BLOCK;
lastBlockWarning = false;
while samplesRead < samplesToRead
    % this isn't really necessary
    assert(ich <= numChannels, 'huh??');

    % read from the file
    % different functions for first vs. subsequent reads
    if firstRead
        blockInfo = PL2ReadFirstDataBlock( filename );
        firstRead = false;
    else
        blockInfo = PL2ReadNextDataBlock( blockInfo );
    end

    % checks

    % block type 2 is analog data.
    %  if there are other block types, we skip and output a warning
    if blockInfo.Type ~= 2
        disp('readPL2Samples: not an analog data block.');
        if lastBlockWarning
            % trim any Nan values that occur in the last block
            naninds = isnan( dataBuffer( 1,: ) );
            % if there are any, they should all be at the end
            if any( naninds ) && sum( naninds > 1 )
                assert( all( diff( find( naninds ) ) == 1 ), ...
                             'readPL2Samples: getting Nans not at end' )
            end

            % trim
            daataBuffer = dataBuffer( :, ~naninds );
            
            break;
        else
            disp('readPL2Samples: haven''t reached last block yet.');
        end
    end
    % should be reading in channels sequentially.
    % if channel # doesn't match expectations, something is very wrong
    assert( blockInfo.Channel == whichChannels( ich ), 'huh??' );
    
    % we expect to always read in NUM_SAMPLES_PLEXON_BLOCK
    %   but this can obv. be shorter at the end of the file
    if blockInfo.AnalogData.NumSamples ~= expectedSamples ...
            && ~lastBlockWarning
        disp( 'readPL2Samples: Warning: did not read in expected number of samples. Hopefully this is the end of the file?' );
        lastBlockWarning = true;
    end

    % store in buffer
    dataBuffer( ich, samplesRead + (1:blockInfo.AnalogData.NumSamples) ) = blockInfo.AnalogData.Values;

    % next channel
    ich = ich+1;
    if ich == numChannels+1
        % next buffer
        ich = 1;
        samplesRead = samplesRead + blockInfo.AnalogData.NumSamples;
    end
end

