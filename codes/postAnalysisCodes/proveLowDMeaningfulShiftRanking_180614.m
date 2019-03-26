%% add your paths here.

% add paths for Feng
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')

addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/postAnalysisCodes')

%% test jPCA code
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools/fromMarksLibraries')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools/CircStat2010d')
%%
buildRuns_20180614

 
%%
loadChoppedCombined_twoLocations


%% make a place to store output videos
outdir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/traj_LowD/meaningful/dim_ranked';
if ~isdir( outdir )
    mkdir( outdir );
end

cd(outdir)

%%
%
% number of trials for each day
numTrialsTot = cellfun( @numel, alf );


% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial


for nday = 1 : numel( alf )
    isCorrectArray{ nday } = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nday }.arrayShapesCorrect);
    trialsToKeep{ nday } = isCorrectArray{ nday };% & UE2.isHoldTrial;

    cueLocs{ nday } = unique(UEs{ nday }.cueLoc);

    for nc = 1 : numel( cueLocs{ nday } )
        trialsByCueLoc{ nday }{nc} = find( trialsToKeep{ nday } & (UEs{ nday }.cueLoc==cueLocs{ nday }(nc)));
        rtsByCueLoc{ nday }{nc} = UEs{ nday }.rt( trialsByCueLoc{ nday }{nc} );    
    end
end

%

% concatenate all the factors

% % this window was used for finding oscillations during arrayDelay in
% %         factors 6 7 8 for session 6
%window = round( [-300  00] / binsize_rescaled );

window = round([40 300]/binsize_rescaled);

%window = round( [-1200 : 500] / binsize );

%whichfieldDimred = 'arrayDim';
%whichfieldDimred = 'arrayOnset';
whichfieldDimred = 'cueOnset';

%whichfieldPlot = 'arrayDim';
%newWindow = round( [-900  650] / binsize_rescaled );



%whichfieldPlot = 'arrayOnset';
%newWindow = round( [-200  1100] / binsize_rescaled );

whichfieldLL = 'cueOnset';
newWindow = round( [40  300] / binsize_rescaled );

%whichfieldPlot = 'cueOnset';
%newWindow = round( [0  600] / binsize_rescaled );

%%
% only do dimred based on day 6
timePoints = window(1):window(2);
numBins = numel( timePoints);
numFactors = size( alf{ nday }(1).rates, 1);
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep(6) ) );
allFactors = zeros( numFactors, numBins * totalTrialsToKeep );


for nday = 6
    ind = 1;
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        allFactors( :, (0:numBins-1) + ind ) = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        ind = ind + numBins;
    end
end

%% get trialsToKeepInds
for nday = 1 : numel( alf )
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
end

%% dimred based on all days
timePoints = window(1):window(2);
numBins = numel( timePoints );
numFactors = size( alf{ nday }(1).rates, 1);
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );
allFactors = zeros( numFactors, numBins * totalTrialsToKeep );

%
ind = 1;
for nday = 1 : numel( alf)
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        allFactors( :, (0:numBins-1) + ind ) = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        ind = ind + numBins;
    end
end 

%% calculate number of neurons for each day
for nday = 1 : numel( alf )
    nNeurons(nday) = size(alf{nday}(1).FR, 1);
end

%% load all weight matrix mapping from factors to LFADS rates
weights = [];
rootPath = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/runs/withExternalInput_20180614/param_SGorjS/all/pbt_run/g061_w19/model_params';
for nday =1: numel(alf)
    dayFile_W = ['/LFADS_' datasets(nday).shortName '_cueOnArrayOnTargetDim_HoldRel.h5_out_fac_linear_W:0'];
    dayFile_b = ['/LFADS_' datasets(nday).shortName '_cueOnArrayOnTargetDim_HoldRel.h5_out_fac_linear_b:0'];
    weights(nday).W = h5read([rootPath], dayFile_W);
    weights(nday).W_T = weights(nday).W';
    weights(nday).b = h5read([rootPath], dayFile_b);
    weights(nday).b_T = weights(nday).b';
end

%% concatenate all dahys' weights together
Wrates_W_cat = [weights.W_T]';
Wrates_b_cat = [weights.b_T]';

%% project from factor space to rate space
allRates = Wrates_W_cat*(allFactors) + Wrates_b_cat;


%% do pca
meanRates = mean( allRates' );

[pca_proj_mat, pc_data, latent] = pca( allRates', 'NumComponents', numFactors);

%% get projected data and spiking data for all trials for cue = 3 and all days (chopped to the desired period)
% set up a cell array for store all post-processed data (spikes, projected data, partial-averaged factors and reconstructed rates)
window = newWindow(1):newWindow(2);
numBins = numel( window );
apd = {};
for nday = 1 : numel( alf )
    itc = 1;
    for itr = 1 : numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        if UEs{nday}.cueLoc( ntr ) == 1
            continue;
        end
        apd{nday}(itc).spikes = alf{nday}(ntr).spikes(:, alf{nday}(ntr).(whichfieldLL) + window);
        apd{nday}(itc).FR = alf{nday}(ntr).FR(:, alf{nday}(ntr).(whichfieldLL) + window);
        frep = alf{ nday }( ntr ).rates( :, alf{ nday }( ntr ).( whichfieldLL ) + window );
        frep = Wrates_W_cat*(frep) + Wrates_b_cat;
        apd{nday}(itc).reconstructedRateSpace = exp(frep); % add FR calculated from LFADS rates to apd, for later use to calculate LL (sanity check)
        frep = frep - repmat( meanRates(:), 1, numBins );

        % project this data
        apd{nday}(itc).rates_lowD = pca_proj_mat' * frep;
        itc = itc + 1;
    end
end


%% iterate through dimensions for averaging
avg_type ='allDays';
%avg_type = 'eachDay';
LL_neurons = {};
for itr = 1:numFactors
    sum_lowD = {};
    for nday = 1:numel(apd)
        sum_lowD{nday} = zeros((numFactors - itr + 1), numBins);
        for tr = 1:numel(apd{nday})
            sum_lowD{nday} = sum_lowD{nday} + apd{nday}(tr).rates_lowD(itr:end,:);
        end
    end
    if strcmp(avg_type, 'allDays')
        sum_avg_dims = zeros((numFactors - itr + 1), numBins);
        for nday = 1:numel(apd)
            sum_avg_dims = sum_avg_dims + sum_lowD{nday};
        end
        avg_dims = sum_avg_dims/sum(cellfun(@numel, apd));
        for nday = 1:numel(apd)
            for tr = 1:numel(apd{nday})
                if itr == 1
                    apd{nday}(tr).avg_lowD = avg_dims;
                else
                    apd{nday}(tr).avg_lowD = cat(1, apd{nday}(tr).rates_lowD(1:(itr - 1), :), avg_dims);
                end
            end
        end
    elseif strcmp(avg_type, 'eachDay')
        for nday = 1:numel(apd)
            avg_dims = sum_lowD{nday}/numel(apd{nday});
            for tr = 1:numel(apd{nday})
                if itr == 1
                    apd{nday}(tr).avg_lowD = avg_dims;
                else
                    apd{nday}(tr).avg_lowD = cat(1, apd{nday}(tr).rates_lowD(1:(itr - 1), :), avg_dims);
                end
            end
        end
    end

    % project the data back to original rate space & pick up the dimensions corresponding to each day
    for nday = 1:numel(apd)
        for tr = 1:numel(apd{nday})
            apd{nday}(tr).rates_trans = pca_proj_mat*inv(pca_proj_mat'*pca_proj_mat)*apd{nday}(tr).avg_lowD + repmat(meanRates(:), 1, numBins);
            %            apd{nday}(tr).FR_trans = 1000/(binsize_rescaled)*exp(weights(nday).W*apd{nday}(tr).factors_trans + weights(nday).b);
            neuronId = (sum(nNeurons(1:(nday-1))) + 1):sum(nNeurons(1:nday));
            apd{nday}(tr).FR_trans = exp(apd{nday}(tr).rates_trans(neuronId,:));
        end
    end
    % For each neuron on each day, calculate a log likehood
    for nday = 1:numel(apd)
        FR_allNeurons = [apd{nday}.FR_trans];
        spiking_allNeurons = [apd{nday}.spikes];
        for n = 1:size(spiking_allNeurons,1)
            FR_thisNeuron = FR_allNeurons(n,:);
            spiking_thisNeuron = spiking_allNeurons(n,:);
            LL_neurons{nday}(n).ll(itr) = sum( log( poisspdf( spiking_thisNeuron(:), FR_thisNeuron(:) ) ) ) / sum(spiking_thisNeuron);
        end
    end
end

%% add in the case when using LFADS firing rate to calculate LL. Also add in the case when using the reconstructed rate space to calculate LL
% these two LL should be equal
for nday = 1:numel(apd)
    %    FR_allNeurons = [apd{nday}.FR];
    FR_allNeurons = [apd{nday}.FR]/250;
    spiking_allNeurons = [apd{nday}.spikes];
    neuronId = (sum(nNeurons(1:(nday-1))) + 1):sum(nNeurons(1:nday)); % for pickup neurons indices for this day
    recons_FR_allNeurons = [apd{nday}.reconstructedRateSpace]; % for calculate LL using factor-constructed rates
    recons_FR_allNeurons = recons_FR_allNeurons(neuronId, :); % for calculate LL using factor-constructed rates
    for n = 1:size(spiking_allNeurons,1)        
        FR_thisNeuron = FR_allNeurons(n,:);
        recons_FR_thisNeuron = recons_FR_allNeurons(n,:);
        spiking_thisNeuron = spiking_allNeurons(n,:);
        LL_neurons{nday}(n).ll(numFactors + 1) = sum( log( poisspdf( spiking_thisNeuron(:), recons_FR_thisNeuron(:) ) ) ) / sum(spiking_thisNeuron);
        LL_neurons{nday}(n).ll(numFactors + 2) = sum( log( poisspdf( spiking_thisNeuron(:), FR_thisNeuron(:) ) ) ) / sum(spiking_thisNeuron);
        % stores the LL that was calculated using LFADS rates
    end
end
%%
LL = zeros(1,32);
count = 0;
for nday = 1:6    
    for n = 1:numel(LL_neurons{nday})
        sum_inf = sum(isinf(LL_neurons{nday}(n).ll));
        if sum_inf >0
            continue
        end
        LL = LL + LL_neurons{nday}(n).ll;
        count = count + 1;
    end
end
LL = LL/count;

%% calculate all-neuron avg LL and stderr
count = 1;
LL_1 = zeros(120, 32);
for nday = 1:6    
    for n = 1:numel(LL_neurons{nday})
        LL_1(count, :) = LL_neurons{nday}(n).ll;
        count = count + 1;        
    end
end

mean_ll = zeros(1,32);
error_ll = zeros(1,32);
for l = 1:32
    mean_ll(l) = mean(LL_1(:,l));
    error_ll(l) = std(LL_1(:,l))/sqrt(120);
end

%% plot for each neuron
count = 0;
x_value = 1:numFactors + 2;
for nday = 1:6
    for n = 1:numel(LL_neurons{nday})
        sum_inf = sum(isinf(LL_neurons{nday}(n).ll));
        if sum_inf >0
            continue
        end
        f1 = figure
        plot(LL_neurons{nday}(n).ll, '-o');
        ylabel('log likelihood');
        xlabel('non-avg dimensions');
        title(['Day ' int2str(nday)]);
        count = count + 1;
        print(f1,['Neuron ' int2str(count)], '-dpng');
        %printpdf(f1,int2str(nIndices(n)) )
        close;
    end
end

%% plot all-neuron avg without stderr
f2 = figure
x_value = 1:numFactors + 2;
scatter(x_value, LL, 'filled');
ylabel('log likelihood');
xlabel('non-avg dimensions');
title('Avg Log Likelihood Across Neurons');
print(f2, 'all_neurons_avg', '-dpng');

%% plot all-neuron avg with error bar
f3 = figure
x_value = 1:numFactors + 2;
errorbar(x_value, mean_ll, error_ll, '-o', 'MarkerSize',5, 'MarkerEdgeColor','red','MarkerFaceColor','red');
ylabel('log likelihood');
xlabel('non-avg dimensions');
title('Avg Log Likelihood Across Neurons');
print(f2, 'all_neurons_avg_withErrorbar', '-dpng');

%% sanity check. Check if reconstructed firing rate after PC transform still look similar to the LFADS output
for nday = 1:6    
    for tr = 1:numel(apd{nday})
        neuronId = (sum(nNeurons(1:(nday-1))) + 1):sum(nNeurons(1:nday));
        tmp = pca_proj_mat*inv(pca_proj_mat'*pca_proj_mat)*apd{nday}(tr).rates_lowD + repmat(meanRates(:), 1, numBins);
        apd{nday}(tr).rates_afterPC_transform = exp(tmp(neuronId,:));
    end
end
%%
f1 = figure
s1 = subplot(1,2,1);
plot(250*apd{1}(1).rates_afterPC_transform(1,:));
xlabel('timePoints')
ylabel('rate (per second)')
title('Reconstructed rates after PC transform')
s2 = subplot(1,2,2);
plot(apd{1}(1).FR(1,:), 'r');
xlabel('timePoints')
ylabel('rate (per second)')
title('LFADS Output rates')
suptitle('Day 1 - trial 1 - neuron 1')
set(f1, 'Position', [131 417 1653 439]);

