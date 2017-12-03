function [spikeCounts] = TrainToCounts(spikeTrain, window)
%TrainToCounts - this function bins the spike train into spike counts. 
%   For all the time points, use either 0 or 1 to indicate if there is
%   spike on that time point
         spikeTrain = spikeTrain*1000;
         
         if abs(max(spikeTrain)-window) < 1e-03
             spikeTrain(end) = spikeTrain(end)-0.1;
         end
         flooredTrain = unique(floor(spikeTrain));
         spikeCounts = zeros(1,window);
         spikeCounts(flooredTrain+1)= spikeCounts(flooredTrain+1)+1;
         dbstop if error

end

