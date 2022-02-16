%This function performs the vanilla NMF decomposition training only a
%subset of the matrix W
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

function [W, H, reconstructError] = nmf_sparse_well_done_semi(X, params, varargin)

numberOfComponentsPerSource = params.numberOfComponentsPerSource;
numberOfSources = sum(numberOfComponentsPerSource);
betaDivergence = params.betaDivergence;
numberOfIterations = params.numberOfIterations;
tolChange = params.tolChange;
tolError = params.tolError;

sqrteps = sqrt(eps);

betaDivergenceAux = beta_loss_to_float(betaDivergence);
[n_samples, n_features] = size(X);

avg1 = sqrt(mean(X(:)) / numberOfSources);
H0 = avg1 * randn(numberOfSources, n_features);
H0 = abs(H0);
avg2 = sqrt(mean(X(:)) / numberOfComponentsPerSource(1));
W0Up = abs(avg2 * randn(n_samples, numberOfComponentsPerSource(1)));
% W0Up = ones(n_samples, numberOfComponentsPerSource(1)) * 1e-3;
W0 = [W0Up params.W0];
W0 = W0 ./ sqrt(sum(W0.^2));
flag = checkMatrices(W0, H0);
if flag
    error('Input matrices have negative elements');
end

if size(W0, 2) ~= numberOfSources || size(H0, 1) ~= numberOfSources
    error('Input matrices dimensions do not match');
end

reconstructError = zeros(numberOfIterations, 1);
reconstructError0 = KL_divergence(X, W0, H0, params.mu);
reconstructError(1) = reconstructError0;
sigma = W0*H0;

for i = 2:numberOfIterations
    if params.verbose
        if ~mod(i, 10)
            fprintf('Iteration %i\n', i);
        end
    end
    
    if betaDivergenceAux == 2
        numW = X*H0';
        W = max(0, W0 .* (numW ./ (W0*(H0*H0') + eps(numW))));
        numH = W'*X;
        H = max(0, H0 .* (numH ./ ((W'*W)*H0 + eps(numH))));
        reconstructError(i) = norm(X - W*H, 'fro') / 2;
        
    elseif betaDivergenceAux == 1
        
        numH = W0.'* (X .* sigma.^(betaDivergenceAux-2));
        denH = W0.'*(sigma.^(betaDivergenceAux - 1)) + params.mu(2);
        denH = denH + eps(numH);
        H = H0 .*(numH./denH);
        
        sigma = W0*H;
        wAux = (W0*W0.');
        sigmaAux = sigma.^(betaDivergenceAux - 2).*X;
        
        numW = ((sigmaAux) + (wAux)*(sigma.^(betaDivergenceAux -1)))*H.';
        denW = (sigma.^(betaDivergenceAux - 1) + wAux*sigmaAux)*H.' + eps(numW);
        
        W0Aux = W0 .* (numW ./ denW);
        W = [W0Aux(:,1:numberOfComponentsPerSource(1)) params.W0];
        W = W ./ sqrt(sum(W.^2));
        sigma = W*H;
        reconstructError(i) = KL_divergence(X, W, H, params.mu);
    end
    
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
