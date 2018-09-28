function [ ] = plotTrialData( Rsucc, run, par, trialFinal, predFlag, outputDir, filePrefix, writeFlag )
    
    close all;
    for i = 1:numel( trialFinal )
        trialNum = trialFinal( i );
        %plotPred = { 'cPec', 'latD','mDel','lBic', 'lTri', 'brac', 'flCR', 'exCR' };
        plotPred = { 'aDel' };
        [ streamPlot, plotTime, dataSize ] = Cherian.plot.genPredAlignment( Rsucc, ...
                                                          run, plotPred, ...
                                                          trialNum, predFlag, writeFlag, outputDir, filePrefix, i );
        %[ posPlot, ~, ~ ] = Cherian.plot.genPredAlignment( Rsucc, run, { 'pos' }, ...
        %                                                  trialNum , predFlag );
        %Cherian.plot.plotNeural( run, par, plotTime, dataSize, trialFinal( i ), writeFlag, outputDir, filePrefix, i );
    end
end
