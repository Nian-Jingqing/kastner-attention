function tSNE_computePlot(timeLag, timeAlign, rtVector, assembled, whichType, rtTagSlow, rtTagFast, savedir, spikeBinMs, f1)

    nTrials = length(assembled);
    if strcmp(whichType, 'real')
        nNeurons = size(assembled(1).spike_smoothed, 1);
    elseif strcmp(whichType, 'lfads')
        nNeurons = size(assembled(1).rates, 1);
    end
    
    
    % set a time lag relative to cueOnset for selecting neuron data at a
    % specific time point for analysis;
    % TimeAlign = nTimesLFADS/2; 
    % the time at cueOnset is 800 ms, so it's the 200th in LFADS rates
    TimeForAnalysis = timeAlign + timeLag;
    % get the specific time point that we want to use to predict rt
    RateMatrix = zeros(nTrials, nNeurons);
    % put all the reaction time into a vector
    % pull out cue location for each trial and put into a vector
    if strcmp(whichType, 'real')
       for i = 1:nTrials
           RateMatrix(i,:) = assembled(i).spike_smoothed(:,TimeForAnalysis);
           % pull out the neural data for all neuron at a specific time in that trial,
           % and put them into the ith row in the rate matrix
       end
    elseif strcmp(whichType, 'lfads')
       for i = 1:nTrials
           RateMatrix(i,:) = assembled(i).normDimAlign(:,TimeForAnalysis);
           % pull out the neural data for all neuron at a specific time in that trial,
           % and put them into the ith row in the rate matrix
       end
    end
    timeLabelOnPlot = timeLag*spikeBinMs;
    
    Rate_tsne = tsne(RateMatrix);
    % rtSlowRates = RateMatrix(rtTagSlow,:);
    % Use rtTag to indice trials with slow rt
    % rtFastRates = RateMatrix(rtTagFast,:);
    % Use rtTag to indice trials with fast rt
    %     Rate_tsne = tsne(RateMatrix);
    % apply tsne to reduce the dimensionality of 106 neurons to 2d
    f1
    if strcmp(whichType, 'real')
        s1 = subplot(1,2,1)
        scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
        hold on
        scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
        ylabel('2nd Dimension');
        xlabel('1st Dimension');
        title(s1, 'tSNE on smoothed spiking');
        legend('show')
    elseif strcmp(whichType, 'lfads')
        s2 = subplot(1,2,2)
        scatter(Rate_tsne(rtTagSlow,1),Rate_tsne(rtTagSlow,2),'g','filled','DisplayName','Slow RT');
        hold on
        scatter(Rate_tsne(rtTagFast,1),Rate_tsne(rtTagFast,2),'r','filled','DisplayName','Fast RT');
        ylabel('2nd Dimension');
        xlabel('1st Dimension');
        title(s2, 'tSNE on lfads rates');
        legend('show')
    end
end