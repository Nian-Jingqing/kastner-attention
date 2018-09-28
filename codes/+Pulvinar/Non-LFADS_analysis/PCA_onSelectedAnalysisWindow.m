function [ coeff, score, latent, latent_tot, cum_variance] = PCA_onSelectedAnalysisWindow( condAvgFiringRate, startTime, endTime, nCond, nTimesPerCond )
%PCA_ONSELECTEDANALYSISWINDOW Summary of this function goes here
%   Detailed explanation goes here

windowedCondStruct(nCond).condAvgFiringRateInWindow = 1;

for cond = 1:nCond
    startTimeForCond = (cond-1)*nTimesPerCond + startTime;
    endTimeForCond = (cond-1)*nTimesPerCond + endTime;
    windowedCondStruct(cond).condAvgFiringRateInWindow = condAvgFiringRate( :, startTimeForCond : endTimeForCond );
end

condAvgFiringRateForAnalysis = [windowedCondStruct.condAvgFiringRateInWindow];
%% actually run PCA
[coeff, score, latent] = pca( condAvgFiringRateForAnalysis' );

%% calculate total variance explained
% normalize the latents to get variance explained by each PC
latent_tot = latent / sum(latent);
% calculate the cumulative variance explained
cum_variance = cumsum( latent_tot );

end

