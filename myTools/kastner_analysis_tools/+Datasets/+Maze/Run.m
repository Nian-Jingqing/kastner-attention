classdef Run < LFADS.Run
    properties(Dependent)
        fileShellScriptRunTensorboard % Location on disk where shell script for tensorboard will be written
    end
   
    methods
        function r = Run(varargin) 
           r@LFADS.Run(varargin{:});
        end

        function f = get.fileShellScriptRunTensorboard(r)
            f = fullfile(r.path, 'lfads_run_tb.sh');
        end

        function runTensorboard(r)
            system( sprintf('sh %s', r.fileShellScriptRunTensorboard) );
        end

        function writeTensorboardScript(r, portNum, tmux_session_name)
        % function writeTensorboardScript(r, portNum, tmux_session_name)
            f = r.fileShellScriptRunTensorboard;
            fid = fopen(f, 'w');
            tb_path = '/usr/local/lib/python2.7/dist-packages/tensorflow/tensorboard/tensorboard.py';
            cmd = sprintf('python %s --logdir %s --port=%i', tb_path, r.pathLFADSOutput, portNum);
            cmd = LFADS.Utils.tmuxify_string(cmd, tmux_session_name);
            fprintf(fid, cmd);
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end
        
        function seq = convertDatasetToSequenceStruct( r, dataset, mode )
        % load data if necessary
            if isempty( dataset.data )
                data = dataset.loadData();
            else
                data = dataset.data;
            end

            % figure out which channels to use
            whichChannels = r.params.whichChannels;
            if isempty( whichChannels )
                % use all the channels
                whichChannels = 1 : size( data.r(1).y, 1 );
            end
            % how much time to get around the alignment point?
            startOffset = r.params.startOffset;
            endOffset = r.params.endOffset;

            % iterate over trials
            for nr = 1:numel( data.r )
                % which alignment point?
                switch lower( r.params.align )
                  case 'flyappears'
                    alignpoint = data.r( nr ).actualFlyAppears;
                  case 'movementonset'
                    alignpoint = data.r( nr ).offlineMoveOnsetTime;
                  case 'gocue'
                    alignpoint = data.r( nr ).actualLandingTime;
                  otherwise
                    disp(sprintf(' Don''t know how to align %s data', alignpoint));
                end
                dataWindow = alignpoint + [startOffset endOffset];
                assert( dataWindow(1)>=1, 'bad align1');
                assert( dataWindow(2)>=1, 'bad align2');
                assert( dataWindow(2)>dataWindow(1), 'bad align3');

                if size(data.r( nr ).y, 2) < dataWindow(2)
                    data.r( nr ).y = int16( [] );
                    data.r( nr ).X = [];
                    data.r( nr ).T = nan;
                else
                    data.r( nr ).y = int16( full( data.r( nr ).y( whichChannels, dataWindow( 1 ) : dataWindow( 2 ) ) ) );
                    data.r( nr ).X = data.r( nr ).X(:, dataWindow( 1 ): ...
                                            dataWindow( 2 ) );
                    data.r( nr ).T = size( data.r( nr ).X, 2 );
                end
                data.r( nr ).y_time = startOffset : endOffset;

                data.r( nr ).binWidthMs = data.r( nr ).params.dtMS;
                data.r( nr ).conditionId = data.r( nr ).conditionCode;
            end

            % if requested, select specific trials
            if ~isempty( r.params.whichTrials )
                ttk = numel( r.params.whichTrials );
                alltrials = numel( data.r );
                disp( sprintf( 'removing %g / %g trials, following r.params.whichTrials', ...
                               alltrials - ttk, alltrials ) );
                data.r = data.r( r.params.whichTrials );
            end

            %remove trials that are too short
            if any(isnan( [ data.r.T ] ))
                disp( sprintf( 'warning: removing %i / %i trials', sum( isnan( [data.r.T] ) ),...
                               numel( data.r ) ) );
                seq = data.r( ~isnan( [ data.r.T ] ) );
            else
                seq = data.r;
            end

            % if requested, shorten the total number of trials
            nTrialsKeep = r.params.nTrialsKeep;
            if nTrialsKeep ~= 0 && numel( seq ) > nTrialsKeep
                seq = seq(1:nTrialsKeep);
            end


        end
        
        function makeLFADSInput(r, regenerate)
            % modified to save two sets of output files one is full data,
            % another one is only held-in trials
            % Generate the LFADS input HD5 files and save them to disk in the pathCommonData folder.
            % If a file already exists, keep the existing file unless
            % regenerate is true. Then symlink the HD5 files used by this
            % run into pathLFADSInput.
            %
            % Args:
            %   regenerate (bool) : Regenerate HD5 files on disk. If false,
            %     the existing files will be left alone.

            if nargin < 2
                regenerate = false;
            end

            seqs = {};
            validInds = {};
            trainInds = {};

            par = r.params;
            if regenerate
                r.deleteSequenceFiles();
            end

            % check which files need to be regenerate
            maskGenerate = false(r.nDatasets, 1);
            fnames = r.lfadsInputFileNames;
            inputInfoNames = r.lfadsInputInfoFileNames;
            if ~regenerate
                for iDS = 1:r.nDatasets
                    fname = fullfile(r.pathCommonData, fnames{iDS});
                    if ~exist(fname, 'file')
                        maskGenerate(iDS) = true;
                    end

                    fname = fullfile(r.pathCommonData, inputInfoNames{iDS});
                    if ~exist(fname, 'file')
                        maskGenerate(iDS) = true;
                    end
                end
            else
                maskGenerate = true(r.nDatasets, 1);
            end

            if any(maskGenerate)
                regenerate = false;
                seqData = r.loadSequenceData(regenerate); % this will set r.sequenceData

                r.assertParamsOkayForSequenceData(seqData);

                % no need to regenerate if alignment and export use the
                % same data, since we would have just regenerated them via
                % loadSequenceData above
                regenerateAlignmentData = regenerate && r.usesDifferentDataForAlignment();
                % ### MRK TODO: for multisession alignment matrices must be calculated on
                % held-in conditions only 
                if r.nDatasets > 1 && r.params.useAlignmentMatrix
                    % generate alignment matrices for stitching run 
                    useAlignMatrices = true;
                    [alignmentMatrices, alignmentBiases] = r.doMultisessionAlignment(regenerateAlignmentData);
                    
                elseif r.version >= 20171107 && r.nDatasets == 1 && r.params.useSingleDatasetAlignmentMatrix
                    % generate alignment matrix for single run (just PCA down to c_factors_dim)
                    useAlignMatrices = true;
                    [alignmentMatrices, alignmentBiases] = r.doMultisessionAlignment(regenerateAlignmentData);
                    
                else
                    % no alignment matrices
                    useAlignMatrices = false;
                end

                % choose validation and training trial indices
                [validIndsCell, trainIndsCell] = deal(cell(r.nDatasets, 1));
                for iDS = 1:r.nDatasets
                    allInds = 1:numel(seqData{iDS});
                    validIndsCell{iDS} = 1 : (r.params.trainToTestRatio+1) : numel(seqData{iDS});
                    trainIndsCell{iDS} = setdiff(allInds, validIndsCell{iDS});
                end
                
                % support old .params.dtMS field
                if isfield(seqData{1}(1), 'params') && isfield(seqData{1}(1).params, 'dtMS')
                    inputBinSizeMs = seqData{1}(1).params.dtMS;
                elseif isfield(seqData{1}(1), 'binWidthMs')
                    inputBinSizeMs = seqData{1}(1).binWidthMs;
                else
                    error('Sequence data lacks binWidthMs field');
                end

                % arguments for the 'seq_to_lfads' call below
                seqToLFADSArgs = {'binSizeMs', par.spikeBinMs,  ...
                    'inputBinSizeMs', inputBinSizeMs, ...
                    'trainInds', trainIndsCell(maskGenerate), 'testInds', validIndsCell(maskGenerate)};

                if useAlignMatrices
                    seqToLFADSArgs{end+1} = 'alignment_matrix_cxf';
                    seqToLFADSArgs{end+1} = alignmentMatrices(maskGenerate);
                    seqToLFADSArgs{end+1} = 'alignment_bias_c';
                    seqToLFADSArgs{end+1} = alignmentBiases(maskGenerate);
                end

                % write the actual lfads input file
                LFADS.Utils.mkdirRecursive(r.pathCommonData);
                
                % save all datasets trials
                LFADS.Interface.seq_to_lfads(seqData(maskGenerate), r.pathCommonData, r.lfadsInputFileNames, ...
                    seqToLFADSArgs{:});
                
                % save only held-in trials to another h5 data file
                if numel(r.params.heldout_trials) > 0
                    idx = 1:numel(seqData);
                    idx = idx(maskGenerate); % only keep Generate indices 
                    heldInfnames = r.lfadsInputFileNames;
                    seqToLFADSArgsHeldin = seqToLFADSArgs;
                    for i = idx
                        heldInfnames{i} = ['heldin' heldInfnames{i}];
                        % remove held-out trials from train/valid indices
                        % remove held_out trials from training sets
                        trials_to_remove = ismember(seqToLFADSArgsHeldin{6}{i}, r.params.heldout_trials{i});
                        seqToLFADSArgsHeldin{6}{i}(trials_to_remove) = [];
                        % remove held_out trials from validation sets
                        trials_to_remove = ismember(seqToLFADSArgsHeldin{8}{i}, r.params.heldout_trials{i});
                        seqToLFADSArgsHeldin{8}{i}(trials_to_remove) = [];
                    end
                    LFADS.Interface.seq_to_lfads(seqData(maskGenerate), r.pathCommonData, heldInfnames, ...
                        seqToLFADSArgsHeldin{:});
                end

                % save input info file for each dataset generated
                inputInfoNames = r.lfadsInputInfoFileNames;
                for iDS = 1:r.nDatasets
                    paramInputDataHash = r.params.generateInputDataHash(); %#ok<*NASGU>
                    if maskGenerate(iDS)
                        trainInds = trainIndsCell{iDS};
                        validInds = validIndsCell{iDS};

                        % save time vectors used in the sequence files to
                        % facilitate fast loading of posterior means sampling
                        seq_timeVector = seqData{iDS}.y_time;
                        seq_binSizeMs = inputBinSizeMs;

                        % save the rebinned spike counts and condition ids
                        % too
                        counts = cat(3, seqData{iDS}.y); % nNeurons x nTime x nChannels
                        if isnumeric(seqData{iDS}(1).conditionId)
                            conditionId = cat(1, seqData{iDS}.conditionId);
                        else
                            conditionId = LFADS.Utils.makecol({seqData{iDS}.conditionId});
                        end
                        
                        % optionally include ground truth
                        extra = {};
                        if isfield(seqData{iDS}, 'y_true')
                            truth = seqData{iDS}.y_true;
                            extra = union(extra, 'truth');
                        end
                        if isfield(seqData{iDS}, 'externalInputs')
                            externalInputs = seqData{iDS}.externalInputs;
                            extra = union(extra, 'externalInputs');
                        end

                        fname = fullfile(r.pathCommonData, inputInfoNames{iDS});
                        save(fname, 'trainInds', 'validInds', 'paramInputDataHash', 'seq_timeVector', 'seq_binSizeMs', 'conditionId', 'counts', extra{:});
                    end
                end
            end

            % check which files need to be symlinked from pathCommonData
            % into pathLFADSInput
            LFADS.Utils.mkdirRecursive(r.pathLFADSInput);
            maskLink = true(r.nDatasets, 1);
            fnames = r.lfadsInputFileNames;
            fnamesInputInfo = r.lfadsInputInfoFileNames;
            for iDS = 1:r.nDatasets
                % make link relative (from link location in
                % runCollection/param_HASH/runName/lfadsInput/)
                if r.version < 20171107
                    % to % runCollection/data_HASH/file.h5
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnames{iDS});
                else
                    % % runCollection/data_HASH/run_name/file.h5
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), r.name, fnames{iDS});
                end

                linkName = fullfile(r.pathLFADSInput, fnames{iDS});
                if ~exist(linkName, 'file') || regenerate
                    LFADS.Utils.makeSymLink(origName, linkName, false);
                end

                % make link relative
                if r.version < 20171107
                    % runCollection/data_HASH/inputInfo.mat
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnamesInputInfo{iDS});
                else
                    % runCollection/data_HASH/run_name/inputInfo.mat
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), r.name, fnamesInputInfo{iDS});
                end
                linkName = fullfile(r.pathLFADSInput, fnamesInputInfo{iDS});
                if ~exist(linkName, 'file') || regenerate
                    LFADS.Utils.makeSymLink(origName, linkName, false);
                end
            end
            
            % DO THIS FOR HELDIN DATA as well
            % check which files need to be symlinked from pathCommonData
            % into pathLFADSInput
            if numel(r.params.heldout_trials) > 0
                fnames = r.lfadsInputFileNames;
                LFADS.Utils.mkdirRecursive(r.pathLFADSInput);
                maskLink = true(r.nDatasets, 1);
                fnames = heldInfnames;
                for iDS = 1:r.nDatasets
                    % make link relative (from link location in
                    % runCollection/param_HASH/runName/lfadsInput/)
                    if r.version < 20171107
                        % to % runCollection/data_HASH/file.h5
                        origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnames{iDS});
                    else
                        % % runCollection/data_HASH/run_name/file.h5
                        origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), r.name, fnames{iDS});
                    end

                    linkName = fullfile(r.pathLFADSInput, fnames{iDS});
                    if ~exist(linkName, 'file') || regenerate
                        LFADS.Utils.makeSymLink(origName, linkName, false);
                    end
% skip the info file for heldin data
%                     % make link relative
%                     if r.version < 20171107
%                         % runCollection/data_HASH/inputInfo.mat
%                         origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnamesInputInfo{iDS});
%                     else
%                         % runCollection/data_HASH/run_name/inputInfo.mat
%                         origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), r.name, fnamesInputInfo{iDS});
%                     end
%                     linkName = fullfile(r.pathLFADSInput, fnamesInputInfo{iDS});
%                     if ~exist(linkName, 'file') || regenerate
%                         LFADS.Utils.makeSymLink(origName, linkName, false);
%                     end
                end
            end
        end
        
    end
end