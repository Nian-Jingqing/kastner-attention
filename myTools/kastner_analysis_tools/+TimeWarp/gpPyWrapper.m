classdef gpPyWrapper <  matlab.mixin.Copyable
% General Purpose Python Function MATLAB Wrapper
    properties
        output_data_path
        output_handler
        py_cmd
        py_dir
        py_file
        py_name
        sys_args
    end
    methods
        % constructor
        function pyWrap = gpPyWrapper( pyName, pyDir, sysArgs, outputHandler )
            assert( iscell( sysArgs ), 'ERROR: sys_args must be passed as cell array.' )
            % python file name
            pyWrap.py_name = pyName;
            % directory housing python file
            pyWrap.py_dir = pyDir;
            % executable filepath for python script
            pyWrap.py_file = fullfile( pyDir, pyName );
            % system arguments to be passed to python function
            pyWrap.sys_args = sysArgs;
            % generate python command
            pyWrap.genPyCmd()
            % generate output path
            pyWrap.genOutPath()
            % output handler function
            pyWrap.output_handler = outputHandler;

        end
        function [] = genOutPath( pyWrap )
            outputDir =  pyWrap.sys_args{ 3 };
            outputName = [ pyWrap.sys_args{ 4 } '.mat' ];
            pyWrap.output_data_path = fullfile( outputDir, outputName );
        end
        
        function [] = genPyCmd( pyWrap )
            pyCmd = 'python';
            for i=1:numel( pyWrap.sys_args ) + 1
                if i == 1
                    pyCmd = [ pyCmd, ' ', pyWrap.py_file ];
                else
                    pyCmd = [ pyCmd, ' ', pyWrap.sys_args{ i - 1 } ];
                end
            end
            % executable python command 
            pyWrap.py_cmd = pyCmd;
        end
        function [] = updateSysArgs( pyWrap, newSysArgs )
            assert( iscell( newSysArgs ), 'ERROR: sys_args must be passed as cell array.' )
            pyWrap.sys_args = newSysArgs;
            % re-generate python command
            pyWrap.genPyCmd()
            % re-generate output path
            pyWrap.genOutPath()
        end
        
        function [] = call( pyWrap )
        % add directory housing python script to system path to allow for python call anywhere in terminal
            if count( py.sys.path, pyWrap.py_dir ) == 0
                insert( py.sys.path,int32(0), pyWrap.py_dir );
            end

            o = system( pyWrap.py_cmd );
            assert( ~o, sprintf( 'ERROR: %s did not run succesfully\n', pyWrap.py_name ) );
        end
        function [ output ] = fetchOutput( pyWrap )
            try
                output = pyWrap.output_handler( pyWrap.output_data_path );
            catch ME
                disp( 'ERROR: Output has not been generated. Run call() first.')
            end
        end
    end
end
