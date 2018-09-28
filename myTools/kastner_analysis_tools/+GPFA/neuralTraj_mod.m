function result = neuralTraj_mod(runIdx, dat, varargin)
% modification by chethan pandarinath
% simple modification which allows specification of training and
% testing data
%
% result = neuralTraj(runIdx, dat, ...)
%
% Prepares data and calls functions for extracting neural trajectories.
%
% INPUTS:
%
% runIdx      - results files will be saved in mat_results/runXXX, where
%               XXX is runIdx
% dat         - structure whose nth entry (corresponding to the nth experimental
%               trial) has fields
%                 trialId -- unique trial identifier
%                 spikes  -- 0/1 matrix of the raw spiking activity across 
%                            all neurons.  Each row corresponds to a neuron.  
%                            Each column corresponds to a 1 msec
%                            timestep.
% trainTrialIdx
% testTrialIdx    - training and testing trial indexes into the dat struct
%
% OUTPUTS:
%
% result      - structure containing all variables saved in mat_results/runXXX/
%               if 'numFolds' is 0.  Else, the structure is empty.
%               
% OPTIONAL ARGUMENTS:
%
% method      - method for extracting neural trajectories
%               'gpfa' (default), 'fa', 'ppca', 'pca'
% binWidth    - spike bin width in msec (default: 20)
% numFolds    - number of cross-validation folds (default: 0)
%               0 indicates no cross-validation, i.e. train on all trials.
% xDim        - state dimensionality (default: 3)
% saveDir     - directory to output results to
%
% @ 2009 Byron Yu         byronyu@stanford.edu
%        John Cunningham  jcunnin@stanford.edu

  method        = 'gpfa'; 
  binWidth      = 20; % in msec
  xDim          = 8;
  saveDir       = pwd;
  trainTrialIdx = []; % added by CP
  testTrialIdx = []; % added by CP
  extraOpts     = GPFA.Util.assignopts(who, varargin);

  fprintf('\n---------------------------------------\n');

  outDir = fullfile(saveDir,'mat_results');
  if ~isdir(outDir)
      fprintf('\ncreating directory %s\n', outDir);
      mkdir(outDir); %'mat_results');
  end
  % Make a directory for this runIdx if it doesn't already exist
  runDir = sprintf('%s/run%03d', outDir,runIdx);
  if isdir(runDir)
    fprintf('Using existing directory %s...\n', runDir);
  else
    fprintf('Making directory %s...\n', runDir);
    mkdir(runDir);
  end

  % Obtain binned spike counts
  seq  = GPFA.Util.getSeq(dat, binWidth, extraOpts{:});
  if isempty(seq)
    fprintf('Error: No valid trials.  Exiting.\n');
    result = [];
    return;
  end

  % Set cross-validation folds 
  N    = length(seq);

  fprintf('\n===== Cross-validation =====\n');

  % Specify filename where results will be saved
  fname = sprintf('%s/%s_xDim%02d', runDir, method, xDim);

  if exist([fname '.mat'], 'file')
      fprintf('%s already exists.  Skipping...\n', fname);
  end
    
  seqTrain      = seq( trainTrialIdx );
  seqTest       = seq( testTrialIdx );
    
    % Remove inactive units based on training set
    hasSpikesBool = (mean([seqTrain.y], 2) ~= 0);
  
    for n = 1:length(seqTrain)
      seqTrain(n).y = seqTrain(n).y(hasSpikesBool,:);
    end
    for n = 1:length(seqTest)
      seqTest(n).y = seqTest(n).y(hasSpikesBool,:);
    end

    % Check if training data covariance is full rank
    yAll = [seqTrain.y];
    yDim  = size(yAll, 1);
    
    if rank(cov(yAll')) < yDim
      fprintf('ERROR: Observation covariance matrix is rank deficient.\n');
      fprintf('Possible causes: repeated units, not enough observations.\n');
      fprintf('Exiting...\n');
      return
    end
    
    fprintf('Number of training trials: %d\n', length(seqTrain));
    fprintf('Number of test trials: %d\n', length(seqTest));
    fprintf('Latent space dimensionality: %d\n', xDim);
    fprintf('Observation dimensionality: %d\n', sum(hasSpikesBool));

    % If doing cross-validation, don't use private noise variance floor.  
    extraOpts = {extraOpts{:}, 'minVarFrac', -Inf};      

    % The following does the heavy lifting.
    if isequal(method, 'gpfa')
      GPFA.Core_GPFA.gpfaEngine(seqTrain, seqTest, fname,... 
      'xDim', xDim, 'binWidth', binWidth, extraOpts{:});
    
    elseif ismember(method, {'fa', 'ppca', 'pca'})
      GPFA.Core_twostage.twoStageEngine(seqTrain, seqTest, fname,... 
      'typ', method, 'xDim', xDim, 'binWidth', binWidth, extraOpts{:});
    end

    if exist([fname '.mat'], 'file')
      save(fname, 'method', 'hasSpikesBool', 'extraOpts', '-append');
    end
  
  result = [];  
  if (nargout == 1) & exist([fname '.mat'], 'file')
    result = load(fname);
  end
