%% condition-average the factors
% windows for each period
periods = {'cueOnset','arrayOnset','arrayDim'};
windows = [[-200; 600] [-100; 500] [-60; 300]];
winBins = round(windows ./binsize_rescaled);
nFactors = size(alf{1}(1).rates, 1);
totBins = sum(diff(winBins)+1);

nConditions = numel( cueLocs{1} );
avgFactors = zeros(nFactors, totBins * nConditions);

%where in the avgFactors matrix are we:
ind = 1;

%iterate over conditions
for ncond = 1:nConditions
    % iterate over task phases
    for nwin = 1:numel( periods )
        % timepoints for this phase
        timepoints = winBins(1, nwin): winBins(2,nwin);
        numBins = numel(timepoints);
        % how many trials did we use for this condition/phase
        totalTrialsThisConditionPhase = 0;

        %iterate over days
        for nd = 1:numel( alf )
            % different trials are needed depending on the phase
            if nwin == 1
                trialsToKeepInds{ nd } = find( UEs{ nd }.cueLoc == cueLocs{ nd }( ncond ) );
            else
                isCorrectArray = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nd }.arrayShapesCorrect);
                trialsToKeepInds{ nd } = find( isCorrectArray & ( UEs{ nd }.cueLoc==cueLocs{ nd }( ncond ) ) );
            end

            for itr = 1:numel( trialsToKeepInds{ nd } )
                ntr = trialsToKeepInds{ nd }( itr );
                totalTrialsThisConditionPhase = totalTrialsThisConditionPhase + 1;
                facData = alf{ nd }( ntr ).rates( :, alf{ nd }( ntr ).( periods{ nwin} ) + timepoints );
                avgFactors( :, (0:numBins-1) + ind ) = avgFactors( :, (0:numBins-1) + ind ) +...
                    facData;
            end

        end

        % normalize by number of trials used
        avgFactors( :, (0:numBins-1) + ind ) = avgFactors( :, (0:numBins-1) + ind) ./ totalTrialsThisConditionPhase;
        ind = ind+numBins;
    end
end


%% do pca
meanFactors = mean( avgFactors' );
avgFactors = avgFactors - repmat( meanFactors(:), 1, totBins * nConditions);

%[pca_proj_mat, pc_data] = pca( avgFactors');%, 'NumComponents', 15);
[pca_proj_mat, pc_data] = pca( avgFactors', 'NumComponents', 12);

