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

function [W, H, reconstructError] = nmf_vanilla_orthogonal_W(X, params, varargin)

numberOfSources = params.numberOfSources;
init = params.init;
betaDivergence = params.betaDivergence;
numberOfIterations = params.numberOfIterations;
tolChange = params.tolChange;
tolError = params.tolError;

sqrteps = sqrt(eps);

betaDivergenceAux = beta_loss_to_float(betaDivergence);
[n_samples, n_features] = size(X);

onesMtx = ones(size(X));
onesMtx2 = ones(numberOfSources);

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

if betaDivergenceAux == 2
   reconstructError0 =  norm(X - W0*H0, 'fro') / 2;
elseif betaDivergenceAux == 1
    reconstructError0 = KL_divergence(X, W0, H0, 'ort', params.mu_orth);
end
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
        
        if ~(params.transpose)
            if params.transform
                numW = (X./(W0*H0 + eps(X)))*H0.' + params.mu_orth(1)*W0;
                denW = (onesMtx*H0.') + eps(numW) + params.mu(1) + params.mu_orth(1)*W0*onesMtx2;

                W = W0 .* (numW ./ denW);
            else
                W = W0;
            end
            numH = W.'*(X./(W*H0 + eps(X))) + params.mu_orth(2)*H0;
            denH = W.'*onesMtx + eps(numH) + params.mu(2) + params.mu_orth(2)*onesMtx2*H0;
            
            H = H0 .* (numH ./ denH);
        else
            if params.transform
                numH = W0.'*(X./(W0*H0 + eps(X)));
                denH = W0.'*onesMtx + eps(numH) + params.mu(2);

                H = H0 .* (numH ./ denH);
            else
                H = H0;
            end
            numW = (X./(W0*H + eps(X)))*H.';
            denW = (onesMtx*H.') + eps(numW) + params.mu(1);
            
            W = W0 .* (numW ./ denW);
        end
       
        reconstructError(i) = KL_divergence(X, W, H, 'ort', params.mu_orth);
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