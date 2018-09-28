function model = fitGLM( data, latentField, spikingField, spikingChannel )
% function model = fitGLM( data, latentField, spikingField, spikingChannel )
% 
%
% CP: script to fit simple GLM with poiss link to single channel
%
% data is a trial-ized struct
%   data(n).( latentField ) = latent data (smooth), dims of numLatents x T
%   data(n).( spikingField ) = spiking data to fit to GLM, dims of numChannels x T
%
% spikingChannel specifies which channel to use for fit
%

L = [data.( latentField )];

Z = [ data.( spikingField) ];

if ~exist('spikingChannel', 'var')
    assert( size(Z, 1) == 1, 'must specify a channel to use if spiking data is multi-dimensional' );
    spikingChannel = 1;
end
Z = Z( spikingChannel, :);

% glmfit expects 
%      predictors (here, latents) to be n_observations x n_predictors
%      responses (here, binned spike count) to be n_observations x 1

model.betas = glmfit( L', Z(:), 'poisson');
