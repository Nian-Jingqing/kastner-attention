classdef RunCollection < PBTExperiment.RunCollection
    % no need to modify anything here, but feel free to add useful methods
    % and properties as useful
    
    methods
        function rc = RunCollection(rootPath, name, datasetCollection, varargin)
            rc@PBTExperiment.RunCollection(rootPath, name, datasetCollection, varargin{:});
        end
    end
end