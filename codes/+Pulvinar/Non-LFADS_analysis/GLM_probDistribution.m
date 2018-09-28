%% build the dataset collection

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

%% Locate and specify the datasets
datasetPath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
    'multi-unit/preAligned/multi-day_CoAoTdHoldRel_JanToApr/withGoodNeurons_HoldRelSepForAO_issueFixed/'];
dc = Pulvinar.DatasetCollection(datasetPath);
dc.name = 'multiDay_CO_AO_TD_HoldRelSepForAO_JanToApr';

% add individual datasets
Pulvinar.Dataset(dc, '170127_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170130_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170201_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170211_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170308_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170311_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170320_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170324_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170327_cueOnArrayOnTargetDim_HoldRel.mat');
Pulvinar.Dataset(dc, '170329_cueOnArrayOnTargetDim_HoldRel_1st.mat');
Pulvinar.Dataset(dc, '170329_cueOnArrayOnTargetDim_HoldRel_2nd.mat');
Pulvinar.Dataset(dc, '170331_cueOnArrayOnTargetDim_HoldRel_1st.mat');
Pulvinar.Dataset(dc, '170331_cueOnArrayOnTargetDim_HoldRel_2nd.mat');
Pulvinar.Dataset(dc, '170407_cueOnArrayOnTargetDim_HoldRel.mat');
% add more datasets here if needed, same code

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Post-running analysis - loading data and the output of LFADS

neuronTotal = zeros(1, 14);
predictableNeuron = zeros(1, 14);
for day = 1:14
    
    %% loading and put real data into Rstruct
    r_real = dc.datasets(day).loadData(); % get the original dataset (for all neurons)
    r_real = R.Rstruct(r_real.R); % put the dataset into R struct class

    %% get dataset info

    nTrials = length(r_real.r); % get trial number
    nTimesRaw = size(r_real.r(1).spikeCounts, 2); % get trial length for raw data, AKA, before re-binned
    nNeurons = size(r_real.r(1).spikeCounts, 1); % get neuron nubmer
    %% smooth the raw spiking
    sigma = 20;
    r_real.smoothFieldInR( 'spikeCounts', 'spike_smoothed', sigma, 1);

    %%
    savedirOne = ['/snel/share/share/derived/kastner/nonLFADS_analysis/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr' ...
        '/GLM/probDistribution/predictedFRVsMean/170127/rebinned20ms_smoothed20ms/'];
    cd(savedirOne);

    %%

    ll = zeros(2, nNeurons);
    % ll_shuffle = zeros(1, nNeurons);
    neuronMeetCriteria = 0;
    for channelToTest = 1:nNeurons
        rcopy = r_real.copy();
        clear legend
        % assume our "factors" are the smoothed neural data of neurons 1-105
        % fit a GLM for neuron 106

        keepChannels = true(nNeurons, 1);
        keepChannels( channelToTest ) = false;
        keepChannels = find( keepChannels );

        for itrial = 1:numel( rcopy.r )
            rcopy.r( itrial ).spike_smoothed = rcopy.r( itrial ).spike_smoothed( keepChannels, : );
            rcopy.r( itrial ).spikeCounts = rcopy.r( itrial ).spikeCounts( channelToTest, : );
        end

        %%
        %  bin the data
        rbinned = rcopy.binData( { 'spikeCounts', 'spike_smoothed' }, [20, 20] );

        %%
        %  fit a GLM model
        numTrials = numel(rbinned);
        trainTrials = true( size( rbinned ) );
        validTrials = 1:5:numTrials;
        trainTrials( validTrials ) = false;
        trainTrials = find( trainTrials );

        model = GLM.fitGLM( rbinned( trainTrials ), 'spike_smoothed', 'spikeCounts', 1 );

        %%
        % test our GLM model
        [dataOut, ll_rawSpiking] = GLM.evalGLM( model, rbinned( validTrials ), 'spike_smoothed', 'spikeCounts', 1, 'firingRate' );

        %%
        % plot the predicted firing rates vs the actual binned spiking activity
        binnedSpikes = [ dataOut.spikeCounts ];
        predictedFR = [ dataOut.firingRate ];





        %##################shuffle###############################

        % shuffle the spiking activity for many time and get many shuffled
        % spiking series

        nShuffles = 100;
        ll_shuffled = zeros(1, nShuffles);
        for s = 1:nShuffles
            shuffle_indices = randperm(length(binnedSpikes));
            spiking_shuffled = binnedSpikes(shuffle_indices);
            ll_shuffled(s) = sum( log( poisspdf( spiking_shuffled(:), predictedFR(:) ) ) ) / sum(spiking_shuffled);
        end

        % calulate if the held-out neuron meet the criteria of p < 0.0001
        mean_shuffledDist = mean(ll_shuffled);
        std_shuffledDist = std(ll_shuffled);
        z_score_heldOut = (ll_rawSpiking - mean_shuffledDist)/std_shuffledDist;
        if z_score_heldOut >= 3.71
            neuronMeetCriteria = neuronMeetCriteria + 1;
        end

    %     % fit probability density function to the shuffled distribution
    %     f1 = figure;
    %     h = histfit(ll_shuffled);
    %     delete(h(1));
    %     %set(gca, 'Legend', 'shuffled spiking');
    %     legend('shuffled spiking')
    %     hold on
    %     
    %     y_range = 0 : 0.01 : max(h(2).YData);
    %     x_value = ll_rawSpiking * ones(size(y_range));
    %     plot(x_value, y_range, 'b', 'DisplayName','true spiking');
    %     legend('show');
    %     title( sprintf( 'Multi-unit %g', channelToTest ) );
    %     ylabel('Distribution');
    %     xlabel('Log likelihood');
    %     print(f1,['Multi-unit ' int2str(channelToTest)], '-dpng');
    %     close;

        %#################end of shuffle############################

        %#################plotting predictedFR vs raw spiking################    
    %     f1 = figure;
    %     scatter( predictedFR, binnedSpikes, 10, 'b','filled' );
    % 
    %     title( sprintf( 'Multi-unit %g', channelToTest ) );
    %     ylabel('actual spiking');
    %     xlabel('predicted FR');
    %     print(f1,['Multi-unit ' int2str(channelToTest)], '-dpng');
    %     close;
        %#################end of plotting predictedFR vs raw spiking################


        %#################plotting predictedFR vs Mean####################
    %     mean_rawSpiking = mean(double(binnedSpikes(:)));
    %     mean_series = mean_rawSpiking + zeros(size(binnedSpikes(:)));
    %     lambda = 0:0.01:1;
    %     ll_inter = zeros(1, length(lambda));
    %     for i = 1:length(lambda)
    %         FR_modified = (predictedFR(:) - mean_rawSpiking) * lambda(i) + mean_series;
    %         ll_inter(i) = sum( log( poisspdf( binnedSpikes(:), FR_modified(:) ) ) ) / sum(binnedSpikes);
    %     end
    %     
    %     f1 = figure;
    %     plot(lambda, ll_inter);
    %     title( sprintf( 'LL vs. lambda for Multi-unit %g', channelToTest ) );
    %     xlabel('lambda');
    %     ylabel('Log likelihood');
    %     print(f1,['Multi-unit ' int2str(channelToTest)], '-dpng');
    %     close;
        %##################end of plotting predictedFR vs Mean#################


    end

    predictableNeuron(day) = neuronMeetCriteria;
    neuronTotal(day) = nNeurons;
end





