function gpfaEngine(seqTrain, seqTest, fname, varargin)
%  
% gpfaEngine(seqTrain, seqTest, fname, ...) 
%
% Extract neural trajectories using GPFA.
%
% INPUTS:
%
% seqTrain      - training data structure, whose nth entry (corresponding to 
%                 the nth experimental trial) has fields
%                   trialId (1 x 1)   -- unique trial identifier
%                   y (# neurons x T) -- neural data
%                   T (1 x 1)         -- number of timesteps
% seqTest       - test data structure (same format as seqTrain)
% fname         - filename of where results are saved
%
% OPTIONAL ARGUMENTS:
%
% xDim          - state dimensionality (default: 3)
% binWidth      - spike bin width in msec (default: 20)
% startTau      - GP timescale initialization in msec (default: 100)
% startEps      - GP noise variance initialization (default: 1e-3)
%
% @ 2009 Byron Yu         byronyu@stanford.edu
%        John Cunningham  jcunnin@stanford.edu

  xDim          = 3;
  binWidth      = 20; % in msec
  startTau      = 100; % in msec
  startEps      = 1e-3;
  extraOpts     = GPFA.Util.assignopts(who, varargin);
    
  % For compute efficiency, train on equal-length segments of trials
  seqTrainCut = GPFA.Util.cutTrials(seqTrain, extraOpts{:});
  if isempty(seqTrainCut)
    fprintf('WARNING: no segments extracted for training.  Defaulting to segLength=Inf.\n');
    seqTrainCut = GPFA.Util.cutTrials(seqTrain, 'segLength', Inf);
  end
  
  % ==================================
  % Initialize state model parameters
  % ==================================
  startParams.covType = 'rbf';
  % GP timescale
  % Assume binWidth is the time step size.
  startParams.gamma = (binWidth / startTau)^2 * ones(1, xDim);
  % GP noise variance
  startParams.eps   = startEps * ones(1, xDim);

  % ========================================
  % Initialize observation model parameters
  % ========================================
  fprintf('Initializing parameters using factor analysis...\n');
  
  yAll             = [seqTrainCut.y];
  [faParams, faLL] = GPFA.Core_twostage.fastfa(yAll, xDim, extraOpts{:});
  
  startParams.d = mean(yAll, 2);
  startParams.C = faParams.L;
  startParams.R = diag(faParams.Ph);

  % Define parameter constraints
  startParams.notes.learnKernelParams = true;
  startParams.notes.learnGPNoise      = false;
  startParams.notes.RforceDiagonal    = true;

  currentParams = startParams;

  % =====================
  % Fit model parameters
  % =====================
  fprintf('\nFitting GPFA model...\n');
  
  [estParams, seqTrainCut, LLcut, iterTime] =... 
    GPFA.Core_GPFA.em(currentParams, seqTrainCut, extraOpts{:});
    
  % Extract neural trajectories for original, unsegmented trials
  % using learned parameters
  [seqTrain, LLtrain] = GPFA.Core_GPFA.exactInferenceWithLL(seqTrain, estParams);

  % ==================================
  % Assess generalization performance
  % ==================================
  if ~isempty(seqTest) % check if there are any test trials
    % Leave-neuron-out prediction on test data 
    if estParams.notes.RforceDiagonal
      seqTest = GPFA.Core_GPFA.cosmoother_gpfa_viaOrth_fast(seqTest, estParams, 1:xDim);
    else
      seqTest = GPFA.Core_GPFA.cosmoother_gpfa_viaOrth(seqTest, estParams, 1:xDim);
    end
    % Compute log-likelihood of test data
    [blah, LLtest] = GPFA.Core_GPFA.exactInferenceWithLL(seqTest, estParams);
  end
  
  % =============
  % Save results
  % =============
  vars  = who;
  fprintf('Saving %s...\n', fname);
  keepvars = vars(~ismember(vars, {'yAll', 'blah'}));

  %% CP: get rid of some variables that are not used from here on out, but take up ungodly amounts of memory
  fields2cut = {'Vsm', 'VsmGP'};
  vars2cut = {'seqTrain', 'seqTrainCut', 'seqTest'};

  for nvar = 1:numel(vars2cut)
      for nf = 1:numel(fields2cut)
          containsVar = false;
          tmp = sprintf('containsVar = isfield(%s, ''%s'');', vars2cut{nvar}, fields2cut{nf});
          eval(tmp);
          if containsVar
              tmp = sprintf('%s = rmfield(%s, ''%s'');', vars2cut{nvar}, vars2cut{nvar}, fields2cut{nf});
              eval(tmp);
          end
      end
  end

  % modifications to save really large files (>2GB)
  for nn=1:numel(keepvars), tmp=whos(keepvars{nn}); tmp2(nn)=tmp.bytes; end
  fsize=sum(tmp2);



  % 1 gb is 2^30
  if log2(fsize) > 31
      save(fname, '-v7.3', keepvars{:});
  else
      save(fname, keepvars{:});
  end
