function [W, H, reconstructError] = nmf_sparse_well_done(X, params, varargin)


numberOfSources = params.numberOfSources;
init = params.init;
betaDivergence = params.betaDivergence;
numberOfIterations = params.numberOfIterations;
tolChange = params.tolChange;
tolError = params.tolError;

sqrteps = sqrt(eps);

betaDivergenceAux = beta_loss_to_float(betaDivergence);
[n_samples, n_features] = size(X);

onesMtx = ones(n_samples, n_samples);
switch init
    case 'random'
        avg = sqrt(mean(X(:)) / numberOfSources);
        H0 = avg * randn(numberOfSources, n_features);
        W0 = avg * randn(n_samples, numberOfSources);
        H0 = abs(H0);
        W0 = abs(W0);
        W0 = W0 ./ sqrt(sum(W0.^2));
    case 'custom'
        W0 = params.W0;
        if params.transform
            H0 = params.H0;
        else
            avg = sqrt(mean(X(:)) / numberOfSources);
            H0 = avg * randn(numberOfSources, n_features);
            H0 = abs(H0);
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
reconstructError0 = beta_divergence(X, W0, H0, betaDivergence, true);
reconstructError(1) = reconstructError0;

sigma = W0*H0;

for i = 2:numberOfIterations
    if params.verbose
        if ~mod(i, 10)
            fprintf('Iteration %i\n', i);
        end
    end
    
    numH = W0.'*(X .* (sigma).^ (betaDivergenceAux-2));
    denH = W0.'*(sigma).^(betaDivergenceAux-1);
    H = H0 .* (numH ./ denH);
    
    sigma = W0*H;
    if params.transform
        numW = (sigma.^(betaDivergenceAux-2).*X)*H.' + W0.*(onesMtx*(W0.*(sigma.^(betaDivergenceAux-2)*H.')));
        denW = sigma.^(betaDivergenceAux-1)*H.' + W0 .* (onesMtx*(W0.*((sigma.^(betaDivergenceAux-2).*X)*H.')));
        W = W0 .* (numW ./ denW);
    else
        W = W0;
    end
    W = W ./ sqrt(sum(W.^2));
    sigma = W*H;
    
    reconstructError(i) = beta_divergence(X, W, H, betaDivergence, false);
    dw = max(max(abs(W - W0) / (sqrteps + max(max(abs(W0))))));
    dh = max(max(abs(H - H0) / (sqrteps + max(max(abs(H0))))));
    delta = max(dw, dh);
    
    if delta <= tolChange
        break;
    elseif (reconstructError0 - reconstructError(i)) <= (tolError*max(1, reconstructError0)) && i > 2
        break;
    elseif i == numberOfIterations
        disp('Max iterations reached. Increase this number for better convergence.')
        break
    end
    
    reconstructError0 = reconstructError(i);
    W0 = W;
    H0 = H;
end