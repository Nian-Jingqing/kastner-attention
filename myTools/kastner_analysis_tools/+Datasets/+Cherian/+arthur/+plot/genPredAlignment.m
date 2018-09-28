function [ h, plotTime, dataSize ] = genPredAlignment( Rsucc, run, predType, trialNum, predFlag, writeFlag,  outputDir, filePrefix, trialId )
    addpath( '/snel/home/lwimala/Projects/MATLAB/dependencies/subaxis')
    addpath( '/snel/home/lwimala/bin/analysis_tools' )
    infield = { 'PCA', 'SmoothedSpikes', 'Factors', 'Rates' };
    
    
    % extract trial start index
    tSidx = Rsucc( trialNum ).tSidx;
    % extract trial end index
    tEidx = Rsucc( trialNum ).tEidx;

    % extract original data sample rate
    origSR = run.Data.orig_sR;
    % extract new sample rate
    newSR = run.Data.new_sR;

    % covert sample index to new bin width
    tSidx = Cherian.utils.convertSampleIndex( tSidx, origSR, newSR );
    tEidx = Cherian.utils.convertSampleIndex( tEidx, origSR, newSR );
    
    % determine plot time
    plotTime = tSidx:tEidx;
    dataSize = numel( plotTime );

    targPosition = Cherian.utils.findTargetPosition( Rsucc, trialNum );

    % define colors for plots    
    ratesColor =  [ 51 200 255 ]./255; % bright blue
                                       %ratesColor = [ 11 255 166 ]./ 255; % bright wintergreen    
    factorsColor = [ 255 155 0 ]./255; % bright orange
                                       %factorsColor = [ 65 173 255 ]./255; % bright orange
    pcaColor = [ 255 11 117 ]./255; % salmon pink
    smoothColor = [ 157 12 186 ]./ 255; % bright wintergreen    
    colors = { factorsColor, factorsColor, ratesColor, ratesColor };
    lineStyle = { ':' , '-', ':', '-' };
    % define opacity for plots
    opacity = { 0.8, 0.8, 0.8, 0.8 };
    
    fontSizeLabels = 12;    
    % for each prediction type ie pos, vel, emg
    for iPred = 1 : numel( predType )
        % initialize figure
        h = Plot.blankFigure();
        figSizeX = 10;
        figSizeY = 10;
        % if there are more than one prediction, use subaxis
        if ~strcmp( predType{ iPred }, 'pos' )
            figSizeY = 5;
            %subaxis( numel( predType ), 1, iPred )
        end % if numel
        figure_prop_name = { 'PaperPositionMode','units','Position' };
        figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };
        set( gcf, figure_prop_name, figure_prop_val )

        step = 20;        
        % plot the true signal
        %Cherian.plot.plotTrue( run, plotTime, predType{ iPred }, ...
        %                      targPosition, dataSize, step  );
        hold on
        %box('off')
        % plot each prediction ie pca, spikes, factors, rates
        % plot only smooth spikes and rates
        predIdx = [ 2 4 ];
        labels = run.Data.Analog(5:end,:);
        disp( predType{ iPred } );
        if predFlag
            diffMean = zeros( size( labels, 1 ), 2 );
            for j = 1 : numel( predIdx )

                i = predIdx( j );
                % extract data
                data = run.Data.( infield{ i } );
                % add bias term
                data = [ data; ones( 1, size( data, 2 ) ) ];
                % extract linear fit matrix
                W = run.W{ i };
                % calculate prediction
                pred = W * data;
                diff{ j } =  labels-pred( 5:end-4, : );
                figure()
                Plot.histf( labels( 1, : ), 40, 'FaceColor', [ 0 0 1 ], 'facealpha', 0.5, 'edgecolor', 'none' )
                hold on
                Plot.histf( labels( 3, : ), 40, 'FaceColor', [ 0 1 0 ], 'facealpha', 0.5, 'edgecolor', 'none' )
                Plot.histf( labels( 5, : ), 40, 'FaceColor', [ 0 1 1 ], 'facealpha', 0.5, 'edgecolor', 'none' )
                Plot.histf( labels( 6, : ), 40, 'FaceColor', [ 1 0 0 ], 'facealpha', 0.5, 'edgecolor', 'none' )
                title( 'Histogram of EMG Signals' )
                ylabel( 'Counts' )
                xlabel( 'EMG Value' )
                box off
                Plot.legalpha( 'Pectoralis', 'Ant. Deltoid','Med. Biceps','Lat. Triceps' )
                legend boxoff
                xlim( [ -50 200 ] )
                filename = '180207-cherian-emg-histogram-cPec-aDel-mBic-lTri-signals';
                filepath = fullfile( outputDir, filename );
                print( filepath,'-dpdf','-r300','-bestfit' )
                print( filepath,'-dpng','-r300' )

                diffMean(:, j ) = mean( diff{j}, 2 )./std( labels, 0, 2 );
                for i=1:size(labels, 1 )
                    negIdx = find( labels( i, : ) < 0 );
                    disp( negIdx' )
                    fprintf( '%4i %0.3f\n', numel( negIdx ), mean( labels( i, negIdx ) ) )
                end
                % plot prediction
                %Cherian.plot.plotPred( run, pred, plotTime, ...
                %                      predType{ iPred }, colors{ i }, opacity{ i }, lineStyle{ i }, step );
            end % for j
            disp( diffMean )
            keyboard
        end
        % plot histograms
        figure()
        Plot.histf( diff{1}(1,:), 40, 'FaceColor', [ 255 155 0 ]./255,'facealpha',0.5, 'EdgeColor', 'none' )
        hold on
        Plot.histf( diff{2}(1,:), 40, 'FaceColor', [ 51 200 255 ]./255,'facealpha',0.5, 'EdgeColor', 'none' )
        title( 'Residuals of Pectoralis Prediction' )
        ylim( [ 0 225000 ] )
        xlim( [ -150 150 ] )
        ylabel( 'Counts' )
        xlabel( 'Residual Value' )
        box off
        Plot.legalpha( 'Smoothed Spikes', 'Rates', 'location', 'northwest' )
        legend boxoff
        hold off
        filename = '180207-cherian-emg-histogram-residuals-cPec';
        filepath = fullfile( outputDir, filename );
        print( filepath,'-dpdf','-r300','-bestfit' )
        print( filepath,'-dpng','-r300' )
        figure()
        Plot.histf( diff{1}(3,:), 40, 'FaceColor', [ 255 155 0 ]./255,'facealpha',0.5, 'EdgeColor', 'none' )
        hold on
        Plot.histf( diff{2}(3,:), 40, 'FaceColor', [ 51 200 255 ]./255,'facealpha',0.5, 'EdgeColor', 'none' )
        title( 'Residuals of Ant. Deltoid Prediction' )
        ylim( [ 0 225000 ] )
        xlim( [ -150 150 ] )
        ylabel( 'Counts' )
        xlabel( 'Residual Value' )
        box off
        Plot.legalpha( 'Smoothed Spikes', 'Rates', 'location', 'northwest' )
        legend boxoff
        hold off
        filename = '180207-cherian-emg-histogram-residuals-aDel';
        filepath = fullfile( outputDir, filename );
        print( filepath,'-dpdf','-r300','-bestfit' )
        print( filepath,'-dpng','-r300' )
        figure()
        Plot.histf( diff{1}(5,:), 40, 'FaceColor', [ 255 155 0 ]./255,'facealpha',0.5, 'EdgeColor', 'none' )
        hold on
        Plot.histf( diff{2}(5,:), 40, 'FaceColor', [ 51 200 255 ]./255,'facealpha',0.5, 'EdgeColor', 'none' )
        title( 'Residuals of Med. Biceps Prediction' )
        ylim( [ 0 225000 ] )
        xlim( [ -150 150 ] )
        ylabel( 'Counts' )
        xlabel( 'Residual Value' )
        box off
        Plot.legalpha( 'Smoothed Spikes', 'Rates', 'location', 'northwest' )
        legend boxoff
        hold off
        filename = '180207-cherian-emg-histogram-residuals-mBic';
        filepath = fullfile( outputDir, filename );
        %print( filepath,'-dpdf','-r300','-bestfit' )
        print( filepath,'-dpng','-r300' )
        keyboard
        if writeFlag
            filename = strcat( filePrefix, '-trial-', num2str( trialId ),'-', predType{ iPred } );
            filepath = fullfile( outputDir, filename );
            print( filepath,'-dpdf','-r300','-bestfit' )
            print( filepath,'-dpng','-r300' )
            close all
        end
    end % for iPred
end

%pcaColor = [ 144 71 71 ]./255; % crimson
%pcaColor = [ 179 135 60 ]./255; % burnt orange
%pcaColor = [ 240 207 39 ]./255; % gold
%pcaColor = [ 238 93 32 ]./255; % dull orange
%smoothColor = [ 115 144 71 ]./ 255; % ugly ass green
%smoothColor = [ 83 229 146 ]./ 255; % wintergreen
%smoothColor = [ 65 179 111 ]./ 255; % dull wintergreen
%smoothColor = [ 36 130 176 ]./ 255; % dull blue ( teal )
