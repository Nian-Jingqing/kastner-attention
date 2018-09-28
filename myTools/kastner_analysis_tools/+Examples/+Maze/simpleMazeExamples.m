%%
dpathBase = '/snel/share/share/data/Shenoy/Maze/';

%%
% dataset info:
% J 27 conditions: 2009-09-12, 2009-09-16
% J 108 conditions: 2009-09-18*, 2009-09-25*, 2009-11-06
% 
% N 27 conditions: 2010-08-12*
% N 108 conditions: 2010-09-10*, 2010-09-23*


%% LOAD THE JENKINS 108 CONDITION DATASET
monkey = 'Jenkins';
dcj = Datasets.Maze.DatasetCollection( fullfile( dpathBase, monkey ) );
dcj.monkey = monkey;
dcj.autoDetectDatasets();

% CP: the best-characterized Maze dataset is Jenkins 2009-09-18
%   that is dc.datasets(3)
dcj.datasets(3).loadData();

%% LOAD THE NITSCHKE 108 CONDITION DATASETS
monkey = 'Nitschke';
dcn = Datasets.Maze.DatasetCollection( fullfile( dpathBase, monkey ) );
dcn.monkey = monkey;
dcn.autoDetectDatasets();

%  dcn.datasets 2 & 3 both have 108 conditions
dcn.datasets( 2 ).loadData();
dcn.datasets( 3 ).loadData();

%% PLOT JENKINS HAND TRAJECTORIES
figure( 1 );
set(gcf,'units','pixels');
set(gcf,'color',[1 1 1]);
set(gcf, 'menubar', 'none');
ah = axes('Position', [0 0 1 1]);
dcj.datasets( 3 ).plotTrajectories();

%% PLOT NITSCHKE HAND TRAJECTORIES
figure( 2 );
set(gcf,'units','pixels');
set(gcf,'color',[1 1 1]);
set(gcf, 'menubar', 'none');
ah = axes('Position', [0 0 1 1]);
dcn.datasets( 2 ).plotTrajectories();


%% THE REST OF THE EXAMPLES FOCUS ON JENKINS DATA

%% align and average the jenkins data
r = dcj.datasets(3).data;
aligned = r.alignAndAverage( 'conditionCode', { 'y', 'X' }, [450 450], 'offlineMoveOnsetTime' );

%% % smooth and bin it
% NOTE: this smoothing is not on the continuous data, so there will be horrendous edge effects. 
whichFields = {'y', 'X'};
binsize = 10;
% smoothing window
smoothSigma = 35;
for iCondition = 1:numel( aligned )
    for iVar = 1:numel( whichFields )
        v = whichFields{ iVar };
        data = aligned( iCondition ).( v );

        % smooth the data
        % function [output, k] = gaussianSmooth(input, sigma,inputdt, holdlength, options)
        data = Utils.gaussianSmooth( data', smoothSigma, 1 )';
        
        currentBins = size( data, 2 );
        binsToKeep = floor( currentBins / binsize ) * binsize;

        newNumBins = binsToKeep / binsize;
        data = reshape(data, [size(data, 1) binsize newNumBins]);
        data = squeeze( sum( data, 2 ) );        

        alignedBinned( iCondition ).( v ) = data;
    end
end

%%
% now plot some neurons
for nn = 1:size( alignedBinned(1).y,1 )
    clf;
    for iCondition = 1:numel(alignedBinned)
        plot( alignedBinned( iCondition ).y(nn, 5 : end-5 ) );
        hold on;
    end
    axis( 'tight' );
    title( nn );
    keyboard
end
