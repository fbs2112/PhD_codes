%This function performs the vanilla NMF decomposition
%Inputs:
%X: input spectrogram signal (it must be real)
%params: struct with the following mandatory fields:
%        numberOfSources: defines the decomposition rank
%        init: 'random' or 'custom'
%        betaDivergence: either 'kullback-leibler' or 'euclidean'
%        numberOfIterations: max number of iterations
%        tolChange: tolerance parameter that defines convergence
%        tolError: error parameter that defines convergence
%
%Outputs:
%W: W matrix
%H: H matrix
%ReconstructError: reconstruction error at each iteration
%
%This code was created by Felipe Barboza da Silva
% Copyright (c) School of Engineering. Macquarie University - Australia.
% All rights reserved. 2020.
% DISCLAIMER:
%     This code is strictly private, confidential and personal to its recipients
%     and should not be copied, distributed or reproduced in whole or in part,
%     nor passed to any third party.
% ** Do not duplicate or distribute without written permission from the owners. **

function [W, H, reconstructError] = nmf_sparse_well_done_v2(X, params, varargin)

numberOfSources = params.numberOfSources;
init = params.init;
numberOfIterations = params.numberOfIterations;
tolChange = params.tolChange;
tolError = params.tolError;

sqrteps = sqrt(eps);

betaDivergenceAux = params.betaDivergenceAux;
[n_samples, n_features] = size(X);

switch init
    case 'random'
        avg = sqrt(mean(X(:)) / numberOfSources);
        H0 = avg * randn(numberOfSources, n_features);
        W0 = avg * randn(n_samples, numberOfSources);
        H0 = abs(H0);
        W0 = abs(W0);
    case 'custom'
        if params.transpose && ~params.transform
            H0 = params.H0;
            avg = sqrt(mean(X(:)) / numberOfSources);
            W0 = avg * randn(n_samples, numberOfSources);
            W0 = abs(W0);
        elseif ~params.transpose && ~params.transform
            W0 = params.W0;
            avg = sqrt(mean(X(:)) / numberOfSources);
            H0 = avg * randn(numberOfSources, n_features);
            H0 = abs(H0);
        else
            H0 = params.H0;
            W0 = params.W0;
        end
        
        flag = checkMatrices(W0, H0);
        if flag
            error('Input matrices have negative elements');
        end
        
        if size(W0, 2) ~= numberOfSources || size(H0, 1) ~= numberOfSources
            error('Input matrices dimensions do not match');
        end
        
    otherwise
        error('Initialisation method not supported');
end

reconstructError = zeros(numberOfIterations, 1);
W0 = W0 ./ sqrt(sum(W0.^2));
reconstructError0 = beta_divergence_v2(X, W0, H0, betaDivergenceAux, params.mu);
reconstructError(1) = reconstructError0;
sigma = W0*H0;

for i = 2:numberOfIterations
    if params.verbose
        if ~mod(i, 10)
            fprintf('Iteration %i\n', i);
        end
    end
    
    numH = W0.'* (X .* sigma.^(betaDivergenceAux-2));
    denH = W0.'*(sigma.^(betaDivergenceAux - 1)) + params.mu(2);
    denH = denH + eps(numH);
    
    H = H0 .*(numH./denH);
    if params.transform
        sigma = W0*H;

        wAux = (W0*W0.');
        sigmaAux = sigma.^(betaDivergenceAux - 2).*X;
        numW = ((sigmaAux) + (wAux)*(sigma.^(betaDivergenceAux -1)))*H.';
        denW = (sigma.^(betaDivergenceAux - 1) + wAux*sigmaAux)*H.';
        denW = denW + eps(numW);
        W = W0 .*(numW./denW);
        W = W ./ sqrt(sum(W.^2));
    else
        W = W0;
    end
    sigma = W*H;
    reconstructError(i) = beta_divergence_v2(X, W, H, betaDivergenceAux, params.mu);
    
    dw = max(max(abs(W - W0) / (sqrteps + max(max(abs(W0))))));
    dh = max(max(abs(H - H0) / (sqrteps + max(max(abs(H0))))));
    delta = max(dw, dh);
    
    if delta <= tolChange
        break;
    elseif (reconstructError0 - reconstructError(i)) <= (tolError*max(1, reconstructError0))
        break;
    elseif i == numberOfIterations
        disp('Max iterations reached. Increase this number for better convergence.')
        break
    end
    
    reconstructError0 = reconstructError(i);
    
    W0 = W;
    H0 = H;
end