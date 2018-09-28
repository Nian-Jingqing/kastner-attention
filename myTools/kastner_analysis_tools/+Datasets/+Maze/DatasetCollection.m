classdef DatasetCollection < LFADS.DatasetCollection
    properties
        monkey = ''
    end
    methods
        function ds = DatasetCollection(path)
            ds = ds@LFADS.DatasetCollection(path);
        end
        
        function autoDetectDatasets(dc)
            dc.clearDatasets;
            
            files = dir(dc.path);
            for iF = 1:numel(files)
                % each dataset should be in an "RC" file

                % check that it starts with RC
                if numel(files(iF).name) < 10, continue, end
                if ~strcmp( files(iF).name(1:2), 'RC' ), continue, end
                if ~strcmp( files(iF).name(end-3:end), '.mat'), continue, end

                ds = Datasets.Maze.Dataset(dc, files(iF).name);
                dc.addDataset(ds);

            end
        end
        
        function filterHasHighSNRChannels(dc)
            dc.loadInfo();
            mask = arrayfun(@(ds) ds.nChannelsHighSNR > 0, dc.datasets);
            dc.filterDatasets(mask);
        end
        
        function filterBestSaveTagEachDate(dc)
            % for days with multiple saveTags - to avoid overlap, we only take
            %    the savetag with the most trials
            allDates = arrayfun(@(ds) ds.datenum, dc.datasets);
            [uniqueDays, ~, udAssignments] = unique(allDates);
            
            maskKeep = false(dc.nDatasets, 1);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets);
            for nd = 1:numel(uniqueDays)
                thisDayDatasetInds = find(udAssignments == nd);
                
                numTrials = nTrials(thisDayDatasetInds);
                [~, whichSetToKeep] = max(numTrials);
                maskKeep(thisDayDatasetInds(whichSetToKeep)) = true;
            end
            
            dc.filterDatasets(maskKeep);
        end
        
        function addDataset(dc, ds)
            assert(isa(ds, 'Datasets.Maze.Dataset'), 'Must be a Datasets.Maze.Dataset instance');
            addDataset@LFADS.DatasetCollection(dc, ds);
        end
        
        function t = getDatasetInfoTable(dc)
            dc.loadInfo();
            rowNames = arrayfun(@(ds) ds.name, dc.datasets, 'UniformOutput', false);
            date = arrayfun(@(ds) ds.datestr, dc.datasets, 'UniformOutput', false);
            saveTags = arrayfun(@(ds) strjoin(ds.saveTags, ','), dc.datasets, 'UniformOutput', false);
            nChannels = arrayfun(@(ds) ds.nChannels, dc.datasets, 'UniformOutput', true);
            nChannelsHighSNR = arrayfun(@(ds) ds.nChannelsHighSNR, dc.datasets, 'UniformOutput', true);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets, 'UniformOutput', true);
            
            t = table(date, saveTags, nTrials, nChannels, nChannelsHighSNR, 'RowNames', rowNames);
        end
    end
end
