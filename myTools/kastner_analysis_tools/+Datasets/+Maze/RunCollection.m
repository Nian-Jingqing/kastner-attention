classdef RunCollection < LFADS.RunCollection
    methods
        function rc = RunCollection(varargin)
            rc@LFADS.RunCollection(varargin{:});
        end
        
        function addRun(rc, r)
            assert(isa(r, 'Datasets.Maze.Run'), 'Must be Datasets.Maze.Run instance');
            addRun@LFADS.RunCollection(rc, r);
        end
    end
end