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

function [W, H, reconstructError] = nmf_vanilla_semi(X, params, varargin)

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
%         HAux = H0(1:numberOfComponentsPerSource(1),:);
        
        numW = (X./(W0*H0 + eps(X)))*H0.';
        denW = (ones(size(W0, 1), size(H0, 2))*H0.') + eps(numW) + params.mu(1);
        
        W0Aux = W0 .* (numW ./ denW);
        
        W = [W0Aux(:,1:numberOfComponentsPerSource(1)) params.W0];
        
        numH = W.'*(X./(W*H0 + eps(X)));
        denH = W.'*ones(size(W, 1), size(H0, 2)) + eps(numH) + params.mu(2);
        
        H = H0 .* (numH ./ denH);
        
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

% W0;