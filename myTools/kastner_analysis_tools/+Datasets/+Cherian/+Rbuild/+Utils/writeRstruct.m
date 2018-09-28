function [ ] = writeRstruct( R, filepath, filename, promptSwitch )
% writes R struct to file

% R : R struct for output
% filepath: directory for output
% filename: name of file for output

% get full file name
    outputName = fullfile( filepath, filename );
    % set flag for prompt
    flag = false;
    % until we get a response
    while ~flag
        if promptSwitch
            prompt = [ '[ ', outputName, ' ] : Is this the correct path for output? [y/n]: '];
            % request response
            response = input( prompt, 's' );
            % lower case
            response = lower( response );
            switch response
                % write file
              case 'y'
                flag = true;
                writeR = true;
                disp('INFO: Writing R Struct ...')
                % dont write file, try again and fix names
              case 'n'
                disp( 'ERROR: Fix filepath and filename and try again.' )
                disp( 'ERROR: No file written .... ' )
                disp( 'ERROR: Exiting ...' )
                flag = true;
                writeR = false;
              otherwise
                % repeat request for response
                disp( "ERROR: Could not understand response. Type 'y' or 'n'.")
            end % switch response
        else
            flag = true;
            writeR = true;
            disp('INFO: Output filename ...')
            disp(outputName)
            disp('INFO: Writing R Struct ...')
        end % if promptSwitch
    end % while ~flag
    
    % if specified write R struct
    if writeR
        save( outputName, 'R', '-v7.3' )
    end % if 
    
end % function
