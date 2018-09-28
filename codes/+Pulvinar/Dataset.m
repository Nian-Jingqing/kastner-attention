classdef Dataset < LFADS.Dataset
    methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
            % you might also wish to set ds.name here,
            % possibly by adding a third argument to the constructor
            % and assigning it to ds.name
        end

        function data = loadData(ds)
            % load this dataset's data file from .path
            data = load(ds.path);
        end

        function loadInfo(ds)
            % Load this Dataset's metadata if not already loaded

            if ds.infoLoaded, return; end

            % modify this to extract the metadata loaded from the data file
            data = ds.loadData();
%             ds.subject = 'TargetDim';
            ds.subject = 'CueOnArrayOnTargetDim_HoldRel';
            ds.saveTags = 1;
            ds.datenum  = 20170608;% the data was from 06/08/2017
            
            
%             ds.nChannels = size(data.R(1).spikeCounts, 1);%find the row number of the first trial, and this is the channel number
%             ds.nTrials = size(data.R);% find the trial number

%             ds.nChannels = size(data.r.r(1).spikes, 1);%find the row number of the first trial, and this is the channel number
%             ds.nTrials = size(data.r.r);% find the trial number
            
            
            ds.infoLoaded = true;
        end

    end
end
