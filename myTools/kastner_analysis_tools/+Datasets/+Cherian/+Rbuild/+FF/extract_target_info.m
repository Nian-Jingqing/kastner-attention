function [ target_struct ] = extract_target_info( targ_table, words, mA, mad )
% INPUTS
% targ_table: Table from dataset containing information about target position and time of appearance
% words: words from dataset corresponding to events occuring in the analog stream
% mA: analog data from madstruct

% convert target table to cell struct
targ_cell = table2cell( targ_table );

% get analog start time
try 
    analog_starttime = mA(1).info.starttime;
catch ME
    analog_starttime = mad.analoginfo( 1 ).starttime;
end


% intialize output struct
target_struct = struct();

% find number of trials in target table
n_trials = size( targ_cell, 1 );

% find number of targets in experiment
n_targets = size( targ_cell, 2 );

% for each trials
for i_trial = 1:n_trials
    % for each target in trial
    for i_targ = 1:n_targets

        % offset target time by analog start time
        targ_time = targ_cell{ i_trial, i_targ }( 1 ) - analog_starttime;

        % get x and y position of target
        x_pos = targ_cell{ i_trial, i_targ }( 2 );
        y_pos = targ_cell{ i_trial, i_targ }( 3 );

        % store x and y position as vector
        targ_pos = [ x_pos, y_pos ];

        target{ i_targ } = { targ_time, targ_pos };

    end

    % we know at least one target appeared every trial
    num_targ_appear = 1;

    % determine if other targets appeared during trial
    for i_targ = 2:size( target, 2 )
        if ~isnan( target{ i_targ }{ 1 } )
            num_targ_appear= num_targ_appear + 1;
        end
    end
    
    % we can use the table for determination of number of appeared targets
    % the number of appeared targets is one less than the number of sucessful targets
    % except for the case of three appeared targets, we then have to determine if the last target was successful or not

    if num_targ_appear == 3
        % get first target time
        targ1_time = target{ 1 }{ 1 };

        % find index of word associated with first target time
        targ1_idx = find( words( :, 1 ) == targ1_time );

        % find trial end by searching for first occurrence of trial start word "18" after trial start
        % since we are searching from targ1_idx to end, we have to account for our index value by adding back the targ1_idx
        % we then subtract 2 to find the index of the trial end
        trial_end_idx = find( words( targ1_idx:end, 2 ) == double( 18 ), 1, 'first' ) + targ1_idx - 2;

        % we then use the index to find the word associated with the end of trial
        word_trial_end = int32( words( trial_end_idx, 2 ) );

        % if word equal to trial success word "32" then add another "target appearance"
        % so that we can maintain the same calculation for the number of succesful targets
        if word_trial_end == 32
            num_targ_appear = num_targ_appear + 1;
        end
    end
    
        target_struct( i_trial ).target = target;
        target_struct( i_trial ).numSuccessfulTargets = num_targ_appear - 1;
end

    
            
        
        
    




