classdef Dataset < LFADS.Dataset
    properties
       exclude_neuron_path = ''; 
    end
   methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
            % you might also wish to set ds.name here,
            % possibly by adding a third argument to the constructor
            % and assigning it to ds.name
        end
        function [] = set_exclude_neuron_path( ds, exclude_neuron_path )
            ds.exclude_neuron_path = exclude_neuron_path;
        end

        function r = loadData( ds, iBlock, datasetType )
        % load this dataset's data file from .path
            data = load( ds.path );
            %load( ds.exclude_neuron_path );
            nBlocks = size( data( 1 ).Experiment_blocks, 1 );
            fprintf( 'INFO: Loading experiment block %i of %i from dataset ...\n', iBlock, nBlocks )
            % extract data from mad struct
            if strcmp( datasetType, 'center-out' )
                [ C, trialStruct, startInds, stopInds, exp_name ] = Datasets.Cherian.Rbuild.CO.loadExpBlock( data, iBlock );
            elseif strcmp( datasetType, 'random-walk' )
                [ C, trialStruct, startInds, stopInds, exp_name ] = Datasets.Cherian.Rbuild.FF.loadExpBlock( data, iBlock );
            else
                disp( 'ERROR: datasetType argument not given or not understood. Please specify if the datasetType is "center-out" or "random-walk".' )
            end
            disp( sprintf( 'INFO: Loaded %s', exp_name ) )
            disp( 'INFO: Smoothing fields ...' )
            % smooth fields
            sigma_neural = 90;
            C.smoothField( 'spikes', 'spikes_smoothed', sigma_neural );
            disp( sprintf( 'INFO: Neural Smoothing Sigma: %i', sigma_neural ) )
            % Gaussian smooth emg data
            C.smoothField( 'emg', 'lp_emg_40Hz', 5 );
            C.smoothField( 'emg', 'lp_emg_10Hz', 17 );
            disp( 'INFO: Generating r struct ...' )
            % create r struct
            r = Datasets.Cherian.cherianData( ...
                C.makeTrialsFromData( startInds, stopInds, trialStruct ) );
        end
        function loadInfo(ds)
            % Load this Dataset's metadata if not already loaded

            if ds.infoLoaded, return; end

            % modify this to extract the metadata loaded from the data file
            %data = ds.loadData();
            %ds.subject = data.subject;
            %ds.saveTags = 1;
            %ds.datenum  = datenum(data.datetime);
            %ds.nChannels = size(data.spikes, 2);
            %ds.nTrials = size(data.spikes, 1);

            ds.infoLoaded = true;
        end

    end
end
