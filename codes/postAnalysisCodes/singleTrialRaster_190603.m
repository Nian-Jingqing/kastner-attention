
%%
% number of trials for each day
numTrialsTot = cellfun( @numel, alf );


% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial

minimalDelay = 800;
for nday = 1 : numel( alf )
    isCorrectArray{ nday } = arrayfun(@(x) strcmp(x, 'HRHR'), UEs{ nday }.arrayShapesCorrect);
    isLongDelay{ nday } = ( [ alf{ nday }.arrayDim ] - [ alf{ nday }.arrayOnset] ) > ( minimalDelay/binsize_rescaled ); % try this later
    isCueLoc3{ nday } = UEs{ nday }.cueLoc == 3;
    trialsToKeep{ nday } = isCorrectArray{ nday } & ( isLongDelay{ nday } )' & isCueLoc3{ nday };% & UE2.isHoldTrial;

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
%window = round( [-300  100] / binsize_rescaled );

%window = round([0 400]/binsize_rescaled);

%window = round( [-1200 : 500] / binsize );

newWindow = round( [200 800] / binsize_rescaled );

whichfieldPlot = 'arrayOnset';
%whichfieldDimred = 'arrayOnset';
%whichfieldPlot = 'cueOnset';

%whichfieldPlot = 'arrayDim';
%newWindow = round( [-900  650] / binsize_rescaled );


%whichfieldPlot = 'arrayOnset';
%newWindow = round( [0  1000] / binsize_rescaled );


%whichfieldPlot = 'cueOnset';
%newWindow = round( [0  600] / binsize_rescaled );


timePoints = newWindow(1):newWindow(2);
numBins = numel( timePoints );
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );

%% setup directory to save
%outdir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/reRun180614_20190603/singleTrials%/200_800/';
%if ~isdir( outdir )
%    mkdir( outdir );
%end

cd(outdir)

%% make single trial plots
% ind = 1;
for nday = 1 : numel( alf)
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );

        f1 = figure;        
        % normalize the LFADS rates for each channel
        norm_rates = normalize(alf{ nday }( ntr ).FR(:, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ), 'centered');
        
        sp(1) = subplot(2, 1, 1);
        imagesc( norm_rates );
        colormap(gca, flipud(gray))
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'200','500','800'});
        title(sp(1), 'LFADS rates');
        ylabel('Multi-units');
        %         set(sp(7), 'FontSize', 7);
        %         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
        %         sp(8) = subplot(7,2,14);
        %         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);


        % plot raw spiking
        sp(2) = subplot(2, 1, 2);
        imagesc( alf{ nday }( ntr ).spikes(:, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ) );
        c = flipud(gray);
        c_1 = [c(1:3,:); c(38:64,:)];
        colormap(gca, c_1)
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'200','500','800'});
        title(sp(2), 'Real Spiking');
        ylabel('Multi-units');
        % %         set(sp(8), 'FontSize', 7);


        %         set(sp(1), 'Position', [0.1300 0.5456 0.7750 0.6])
        %         set(sp(2), 'Position', [0.1300 0.1539 0.7750 0.3119]);
        suptitle(['Day ' int2str(nday) ' - trial ' int2str(ntr)]);
        %         set(f1, 'Position', [279 53 648 913]);
        %         set(f1, 'Position', [375 67 1079 899]);
        print(f1, ['Day ' int2str(nday) ' - trial ' int2str(ntr)], '-dpng');
        close;
    end
end 

