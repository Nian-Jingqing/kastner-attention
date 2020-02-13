%%
outdir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/tc_newSP_PBT_191120/singleTrialPlots/tmp/';
if ~isdir( outdir )
    mkdir( outdir );
end
%% calculate cond-avg matrix
events = {'barOn', 'cueOn', 'targetOn'};
conditions = struct;

% define all conditions to analyze
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};

%%
nBar = 1;
nCue = 1;

cond_name = [conditions.cueOn{nCue}, '_', conditions.barOn{nBar}];
%%

% number of trials for each day
numTrialsTot = cellfun( @numel, alf );


% %  trials we want have the UE2.arrayShapesCorrect string 'HRHR'
% %  they must also be hold trials, i.e. UE2.isHoldTrial

for nday = 1 : numel( alf )
    trialsToKeep{ nday } = (UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue);
end

%



newWindow = round( [100 700] / binsize_rescaled );
whichfieldPlot = 'cueOnset';

timePoints = newWindow(1):newWindow(2);
totalTrialsToKeep = sum( cellfun( @sum, trialsToKeep ) );

%% setup directory to save
%outdir = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/reRun180614_20190603/singleTrials%/200_800/';
%if ~isdir( outdir )
%    mkdir( outdir );
%end

cd(outdir)

%% make single trial plots
% ind = 1;
%for nday = 1 : numel( alf)
for nday = 5
    trialsToKeepInds{ nday } = find( trialsToKeep{ nday } );
    for itr = 1:numel( trialsToKeepInds{ nday } )
        ntr = trialsToKeepInds{ nday }( itr );

        f1 = figure;        
        % normalize the LFADS rates for each channel
        norm_rates = normalize(alf{ nday }( ntr ).rates(:, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ), 'centered');
        
        sp(1) = subplot(6, 1, 1);
        imagesc( norm_rates );
        %colormap(gca, flipud(gray))
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'Cue+100','400','700'});
        title(sp(1), 'LFADS rates');
        ylabel('Multi-units');
        %         set(sp(7), 'FontSize', 7);
        %         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
        %         sp(8) = subplot(7,2,14);
        %         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);


        % plot raw spiking
        sp(2) = subplot(6, 1, 2);
        imagesc( alf{ nday }( ntr ).spikes(:, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ) );
        %c = flipud(gray);
        %c_1 = [c(1:3,:); c(38:64,:)];
        %colormap(gca, c_1)
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'Cue+100','400','700'});
        title(sp(2), 'Real Spiking');
        ylabel('Multi-units');
        % %         set(sp(8), 'FontSize', 7);


        sp(3) = subplot(6,1,3);
        plot( alf{ nday }( ntr ).spikes(3, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ), 'r' );
        hold on
        plot(norm_rates(3,:), 'b');
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'Cue+100','400','700'});
        title('neuron 3');
        axis tight

        sp(4) = subplot(6,1,4);
        plot( alf{ nday }( ntr ).spikes(4, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ), 'r' );
        hold on
        plot(norm_rates(4,:), 'b');
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'Cue+100','400','700'});
        title('neuron 4');
        axis tight

        sp(5) = subplot(6,1,5);
        plot( alf{ nday }( ntr ).spikes(5, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ), 'r');
        hold on
        plot(norm_rates(5,:), 'b');
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'Cue+100','400','700'});
        title('neuron 5');
        axis tight

        sp(6) = subplot(6,1,6);
        plot( alf{ nday }( ntr ).spikes(6, alf{ nday }( ntr ).( whichfieldPlot ) + timePoints ), 'r');
        hold on
        plot(norm_rates(6,:), 'b');
        set(gca,'XTick',[1 0.5*numel(timePoints) numel(timePoints)]);
        set(gca,'XTickLabels',{'Cue+100','400','700'});
        title('neuron 6');
        axis tight

        set(gcf, 'Position', [40 53 908 1033])

        %         set(sp(1), 'Position', [0.1300 0.5456 0.7750 0.6])
        %         set(sp(2), 'Position', [0.1300 0.1539 0.7750 0.3119]);
        suptitle(['Day ' dataset(nday).date ' - trial ' int2str(ntr)]);
        %         set(f1, 'Position', [279 53 648 913]);
        %         set(f1, 'Position', [375 67 1079 899]);
        print(f1, ['Day ' dataset(nday).date ' - trial ' int2str(ntr)], '-dpng');
        close;
    end
end 