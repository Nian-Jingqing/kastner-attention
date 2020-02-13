function singleTrial_delay_lowD(alf, UE, binsize_rescaled, dataset)
% this script performs PCA on single trial delay activity and do more advanced analysis based on the lowD
% define all conditions to analyze
%%
barID = [1 2];
cueID = [1 2 3 4 5 6 7 8];
%cueID = 1;
cmap = lines();
clear keepThisTrial trialsToKeepInds trialsByCueLoc allFactors meanFactors
for nday = 1 : numel( alf )
    keepThisTrial{nday} = (ismember(UE{nday}.barType, barID)) & (ismember(UE{nday}.cueType, cueID));
    %keepThisTrial{nday} = (UE{nday}.barType == barID) & (UE{nday}.cueType == cueID);
    trialsToKeepInds{nday} = find(keepThisTrial{nday});
    trialsByCueLoc{nday} = zeros(1, numel(trialsToKeepInds{nday}));
    for itr = 1:numel(trialsToKeepInds{nday})
        ntr = trialsToKeepInds{nday}(itr);
        trialsByCueLoc{nday}(itr) = UE{nday}.cueType(ntr);
    end    
end
%

window = round([-400 0]/binsize_rescaled);
whichfieldDimred = 'targetStart';
whichfieldPlot = 'targetStart';
newWindow = round( [-100  0] / binsize_rescaled );

%window = round([400 700]/binsize_rescaled);
%whichfieldDimred = 'cueOnset';
%whichfieldPlot = 'cueOnset';
%newWindow = round( [300 700] / binsize_rescaled );

timePoints = window(1):window(2);
numBins = numel( timePoints );
numFactors = size( alf{ nday }(1).factors, 1);
totalTrialsToKeep = sum( cellfun( @sum, keepThisTrial ) );
allFactors = zeros( numFactors, numBins * totalTrialsToKeep );

%
ind = 1;
for nday = 1 : numel( alf)
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        allFactors( :, (0:numBins-1) + ind ) = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldDimred ) + timePoints );
        ind = ind + numBins;
    end
end

% perform PCA
meanFactors = mean( allFactors' );
[pca_proj_mat, pc_data, latent, tsquared, explained] = pca( allFactors', 'NumComponents', 10);

%% color coded by conditions
clf
timePoints_new = newWindow(1):newWindow(2);
numBins_new = numel(timePoints_new);
for nday = 1:6
    subplot(2,3,nday)
    day_title = ['Day ' dataset(nday).date];
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        cond = trialsByCueLoc{nday}(itr);

        frep = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints_new );
        % mean center
        frep = frep - repmat( meanFactors(:), 1, numBins_new );

        % project this data
        dim_reduced_data = pca_proj_mat' * frep;
        plot3(dim_reduced_data(1, :), dim_reduced_data(2, :), dim_reduced_data(3, :), 'Color', 0.7*cmap(cond,:), 'LineWidth', 1);
        %plot(dim_reduced_data(2, :), dim_reduced_data(3, :), 'Color', 'k', 'LineWidth', 1);
        hold on
        plot3(dim_reduced_data(1, 1), dim_reduced_data(2, 1), dim_reduced_data(3, 1), 'o', 'MarkerFaceColor', 0.7*cmap(cond,:), 'MarkerEdgeColor', 0.7*cmap(cond,:))
        %plot(dim_reduced_data(2, end), dim_reduced_data(3, end), 'o', 'MarkerFaceColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:), 'MarkerEdgeColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:))
    end
    %set(gca, 'view', [-40.7000, -22.8000]);
    %axis(gca, [ -5 -2 -3 0  -2 2]) % cond1
    axis(gca, [ -5 5 -4 2  -4 4])  % 4 conds
    %axis(gca, [ -4 2 -4 2  -4 4])  % for smoothed factors
    %axis(gca, [ 2 5 -1 2  -1 1])  % for filtered factors
    %axis(gca, [2.3 4.5 -1 1.6]) % for factors, 4 conds, 2D
    title(day_title);
end


%% color coded hit/miss
figure
paxis = [ -4 4 -3 2  -3 4];
timePoints_new = newWindow(1):newWindow(2);
numBins_new = numel(timePoints_new);
for nday = 1:6
    subplot(2,3,nday)
    axis(gca, paxis)
    day_title = ['Day ' dataset(nday).date];
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );

        if alf{nday}(ntr).isMiss == 1 && UE{nday}.isValidTarget(ntr)
            colorHM = 2; % miss, red
        elseif alf{nday}(ntr).isMiss == 0 && UE{nday}.isValidTarget(ntr)
            colorHM = 1; % hit, blue
        else
            continue;
        end
        
        cond = trialsByCueLoc{nday}(itr);

        frep = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints_new );
        % mean center
        frep = frep - repmat( meanFactors(:), 1, numBins_new );

        % project this data
        dim_reduced_data = pca_proj_mat' * frep;


        h = Plot.patchline(dim_reduced_data(1, :)', dim_reduced_data(2, :)', dim_reduced_data(3, :)', 'edgecolor', 0.85*cmap(colorHM, :));

        set(h, 'facealpha', 0.1, 'edgealpha', 0.1);
        %plot(dim_reduced_data(1, 1:end-5), dim_reduced_data(1, 6:end), 'Color', 'k', 'LineWidth', 1);
        hold on
        plot3(dim_reduced_data(1, end), dim_reduced_data(2, end), dim_reduced_data(3, end), 'o', 'MarkerFaceColor', 0.85*cmap(colorHM, :), 'MarkerEdgeColor', 0.85*cmap(colorHM,:))

        
        %plot3(dim_reduced_data(1, :), dim_reduced_data(2, :), dim_reduced_data(3, :), 'Color', 'k', 'LineWidth', 0.5);
        %plot(dim_reduced_data(1, 1:end-5), dim_reduced_data(1, 6:end), 'Color', 'k', 'LineWidth', 1);
        %hold on
        %plot3(dim_reduced_data(1, end), dim_reduced_data(2, end), dim_reduced_data(3, end), 'o', 'MarkerFaceColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:), 'MarkerEdgeColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:))
        %plot(dim_reduced_data(1, end-5), dim_reduced_data(1, end), 'o', 'MarkerFaceColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:), 'MarkerEdgeColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:))
    end
    %set(gca, 'view', [-40.7000, -22.8000]);
    %axis(gca, [ -5 -2 -3 0  -2 2]) % cond1
    axis(gca, [ -4 4 -3 2  -3 4])  % 4 conds
    %axis(gca, [ -4 2 -4 2  -4 4])  % for smoothed factors
    %axis(gca, [ 2 5 -1 2  -1 1])  % for filtered factors
    %axis(gca, [2.3 4.5 -1 1.6]) % for factors, 4 conds, 2D
    title(day_title);
end
set(gcf, 'Position', [44 15 1874 1071])
suptitle('Red - Miss, Blue - Hit')

%% color coded by reaction time
figure
paxis = [ -4 4 -3 2  -3 4];
%clf
timePoints_new = newWindow(1):newWindow(2);
numBins_new = numel(timePoints_new);
for nday = 1:6
    subplot(2,3,nday)
    axis(gca, paxis)
    day_title = ['Day ' dataset(nday).date];
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        
        if alf{nday}(ntr).rt < median(UE{nday}.trialEnd - UE{nday}.targetOn) && ~alf{nday}(ntr).isEarlyError
            colorRT = 2; % fast trials, red
        elseif alf{nday}(ntr).rt >= median(UE{nday}.trialEnd - UE{nday}.targetOn) && alf{nday}(ntr).isMiss == 0
            colorRT = 1; % slow trials, blue
        else
            continue;
        end
        
        cond = trialsByCueLoc{nday}(itr);

        frep = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints_new );
        % mean center
        frep = frep - repmat( meanFactors(:), 1, numBins_new );

        % project this data
        dim_reduced_data = pca_proj_mat' * frep;
        %keyboard
        h = Plot.patchline(dim_reduced_data(7, :)', dim_reduced_data(8, :)', dim_reduced_data(9, :)', 'edgecolor', 0.85*cmap(colorRT, :));

        set(h, 'facealpha', 0.1, 'edgealpha', 0.1);
        %plot(dim_reduced_data(1, 1:end-5), dim_reduced_data(1, 6:end), 'Color', 'k', 'LineWidth', 1);
        hold on
        plot3(dim_reduced_data(7, end), dim_reduced_data(8, end), dim_reduced_data(9, end), 'o', 'MarkerFaceColor', 0.85*cmap(colorRT, :), 'MarkerEdgeColor', 0.85*cmap(colorRT,:))
        %plot(dim_reduced_data(1, end-5), dim_reduced_data(1, end), 'o', 'MarkerFaceColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:), 'MarkerEdgeColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:))
    end
    %set(gca, 'view', [-40.7000, -22.8000]);
    %axis(gca, [ -5 -2 -3 0  -2 2]) % cond1
    %axis(gca, [ -4 4 -3 2  -3 4])  % 4 conds
    %axis(gca, [ -4 2 -4 2  -4 4])  % for smoothed factors
    %axis(gca, [ 2 5 -1 2  -1 1])  % for filtered factors
    %axis(gca, [2.3 4.5 -1 1.6]) % for factors, 4 conds, 2D
    title(day_title);
end
set(gcf, 'Position', [44 15 1874 1071])
suptitle('Red - Fast, Blue - Slow')

%% plot single PCs
nday = 1;
for nPC = 1:10
    subplot(5,2, nPC)
    icount = 0;
    sum_factors = zeros(1,41);
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        cond = trialsByCueLoc{nday}(itr);

        frep = alf{ nday }( ntr ).factors_filtered( :, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints_new );
        % mean center
        frep = frep - repmat( meanFactors(:), 1, numBins_new );

        % project this data
        dim_reduced_data = pca_proj_mat' * frep;
        plot(dim_reduced_data(nPC, :), 'Color', 0.85*cmap(cond, :), 'LineWidth', 1);
        %plot(frep(nPC, :), 'Color', 0.85*cmap(cond, :), 'LineWidth', 1);
        hold on
        sum_factors = sum_factors+frep(nPC,:);
        icount = icount + 1;
    end
    sum_factors = sum_factors/icount;
    %plot(sum_factors, 'k', 'LineWidth', 3);
    title(['PC' int2str(nPC)])
    axis tight
    %ylim([-0.8 4])
end
suptitle(['Day ' dataset(nday).date])
    

%%

figure
xs = 1:1:2000;
a = sin(2*pi*0.006*xs);
b = cos(2*pi*0.006*xs);
hold on
plot(a, b, 'k')
hold on
for nday = 1:12
    %subplot(2,3,nday)
    %plot(a, 'k')
    day_title = ['Day ' dataset(nday).date];
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        cond = trialsByCueLoc{nday}(itr);
        %x = alf{nday}(ntr).targetStart - alf{nday}(ntr).cueOnset;
        x = UE{nday}.targetOn(ntr) - UE{nday}.cueOn(ntr);
        y_a = sin(2*pi*0.006*x);
        y_b = cos(2*pi*0.006*x);
        plot(y_a, y_b, 'o', 'MarkerFaceColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:), 'MarkerEdgeColor', 0.85*cmap(alf{nday}(ntr).isHit+1,:))
        hold on
    end
    title('6Hz');
end

%%

window_tsne = round([60 100]/binsize_rescaled);
whichfieldTsne = 'targetStart';

timePoints_tsne = window_tsne(1):window_tsne(2);
numBins_tsne = numel( timePoints_tsne );
numFactors = size( alf{ nday }(1).factors, 1);
totalTrialsToKeep = sum( cellfun( @sum, keepThisTrial ) );
%allFactors_tsne = zeros( numFactors, numBins_tsne * totalTrialsToKeep );
allFactors_tsne = zeros( 3, numBins_tsne * totalTrialsToKeep );
colorCode = zeros(1, numBins_tsne * totalTrialsToKeep);

%%

ind = 1;
for nday = 1 : numel( alf)
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );
        tmp = alf{ nday }( ntr ).factors( :, alf{ nday }( ntr ).( whichfieldTsne ) + timePoints_tsne );
        tmp = tmp - repmat( meanFactors(:), 1, numBins_tsne );
        allFactors_tsne( :, (0:numBins_tsne-1) + ind ) = pca_proj_mat' * tmp;
        %allFactors_tsne( :, (0:numBins_tsne-1) + ind ) = alf{ nday }( ntr ).factors_smoothed( :, alf{ nday }( ntr ).( whichfieldTsne ) + timePoints_tsne );
        if alf{nday}(ntr).isHit == 1
            colorCode((0:numBins_tsne-1) + ind) = 4;
        else
            colorCode((0:numBins_tsne-1) + ind) = 10;
        end
        %colorCode((0:numBins_tsne-1) + ind) = alf{nday}(ntr).isHit+1;
        ind = ind + numBins_tsne;
    end
end

%%
figure
Y = tsne(allFactors_tsne');

%%
clf
scatter(Y(:,1), Y(:,2), 5, colorCode)

%%
for i = 1:numel(colorCode)
    scatter(Y(i,1), Y(i,2), 'MarkerEdgeColor', 0.85*cmap(colorCode(i), :), 'MarkerFaceColor', 0.85*cmap(colorCode(i), :))
    hold on
end

