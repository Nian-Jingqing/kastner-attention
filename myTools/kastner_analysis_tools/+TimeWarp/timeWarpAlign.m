function [ data_aligned, shifts ] = timeWarpAlign( dataDir, data )
%keyboard
%% generate input and output directories
inputDataDir = fullfile( dataDir, 'input_data' );
outputDataDir = fullfile( dataDir, 'output_data' );

% if folder does not exist, make folder
disp( 'INFO: Generating directories for input and output data ...')

% generate input data directory
[ status, msg, ~ ] = mkdir( inputDataDir );
if ~isempty( msg )
    msg = [ 'INFO: Input' msg ];
else
    msg = 'INFO: Input Directory generated.';
end
disp( msg )

% generate output data directory
[ status, msg, ~ ] = mkdir( outputDataDir );
if ~isempty( msg )
    msg = [ 'INFO: Output' msg ];

else
    msg = 'INFO: Output Directory generated.';
end
disp( msg )

%% generate data files
inName = '170308_raw_spiking';
outName = '170308_shiftOnly';
overwriteFlag = logical( 1 );
TimeWarp.genTimeWarpDataFile( data, inputDataDir, inName, overwriteFlag );
sysArgs = { inputDataDir, inName, outputDataDir, outName };

%% get the directory that stores the python function
pyName = 'shiftOnly_warp.py';
pyDir = '/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools/+TimeWarp/affinewarp/';
outputHandler = @TimeWarp.timeWarpOutHandler;

%% run time warpping function
disp( 'Running time warpping for:' )
disp( sysArgs{ 2 } )
timeWarping = TimeWarp.gpPyWrapper( pyName, pyDir, sysArgs, outputHandler );
%keyboard
%%
if ~overwriteFlag
    try
        outputTW = timeWarping.fetchOutput();
    catch ME
        % *Run time warping*
        timeWarping.call()
        outputTW = timeWarping.fetchOutput();
    end
else
    % *Run time warping*
    timeWarping.call()
    outputTW = timeWarping.fetchOutput();
end
data_aligned = outputTW.alignedData;
shifts = outputTW.shifts;
end