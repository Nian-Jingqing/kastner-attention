function [data ll l0] = evalGLM( model, data, latentField, spikingField, spikingChannel, outputField )
% function [data ll l0] = evalGLM( model, data, latentField, spikingField, spikingChannel, outputField )
% 
%
% CP: script to fit simple GLM with poiss link to single channel
%
% model is the output of +GLM.fitGLM
%
% data is a trial-ized struct
%   data(n).( latentField ) = latent data (smooth), dims of numLatents x T
%   data(n).( spikingField ) = spiking data to fit to GLM, dims of numChannels x T
%
% spikingChannel specifies which channel to use for fit
%

% iterate over trials
for itrial = 1:numel( data )
    % pull out latent variables
    Ltrial = data( itrial ).( latentField );
    % use the function glmval to get a predicted firing rate given the latent variables
    %    and the model
    poisspred = glmval( model.betas, Ltrial', 'log');
    data( itrial ).( outputField ) = poisspred(:)';
end

% get firing rates for all trials
poisspred = [ data.( outputField ) ];

Z = [ data.( spikingField) ];
if ~exist('spikingChannel', 'var')
    assert( size(Z, 1) == 1, 'must specify a channel to use if spiking data is multi-dimensional' );
    spikingChannel = 1;
end
Z = Z( spikingChannel, :);

% evaluate the log likelihood of the observed spiking activity given the predicted firing rate
ll = sum( log( poisspdf( Z(:), poisspred(:) ) ) ) / sum(Z);
% evaluate the log likelihood of the observed spiking activity assuming a constant mean FT
l0 = sum( log( poisspdf( Z(:), mean(double(Z(:))) + zeros(size(Z(:)))) ) ) / sum(Z);
