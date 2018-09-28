%% set your paths
addpath('/snel/home/fzhu23/bin/analysis_tools');

%% get file directory and file names
ddir = '/snel/share/share/data/Kastner_Pulvinar/singleSession/continuous/data_raw/M20170608_PUL_all-g2-g3-g4-evokedSpiking-v8';
fnames = fullfile( ddir, '*.mat' );
theFiles = dir(fnames);

%% sort theFIles in the sequence based on modification time
cells = struct2cell(theFiles); % put the struct array into cell array
sortVals = cells(3,:)'; % take out the modification date from theFiles
ix = (1:length(sortVals))'; % initialize a vector for indexing the file sequence
time = datetime(sortVals); % convert modification time from cell to datetime for sorting
timeTable = table(time,ix); % put time and index in a table
sorted = sortrows(timeTable,'time'); % sort index by modification time 
theFiles = theFiles(sorted.ix); % use the index to sort theFiles

