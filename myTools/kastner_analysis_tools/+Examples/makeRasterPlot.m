%% this is some stupid example code to demonstrate how to make raster plots

%% make some fake data
seq = Examples.makeFakeData();

allConditions = [ seq.conditionID ];
conditions = unique( allConditions );
numPanels = ceil( sqrt( numel( conditions ) ) );

binsize = 0.01;

neuronsToPlot = 1:2;

%iterate over some of the neurons
for nn =  neuronsToPlot
    figure( nn ); clf;

    % separate trials out by condition
    for nc = 1:numel( conditions )
        trials = find( allConditions == conditions( nc ) );

        spks = {};
        % get the spiking data for this neuron, all individual trials
        for nt = 1:numel( trials )
            spkTimes{ nt } = find( seq( trials( nt ) ).y( nn, : ) );
            % the spikes are binned, so scale times appropriately
            spkTimes{ nt } = spkTimes{ nt } * binsize;
        end

        subplot( numPanels, numPanels, nc );
        Plot.tickRaster( spkTimes );

        axis('tight');
        xlim( [ 0 size( seq(1).y, 2) ] * binsize );
    end

end
