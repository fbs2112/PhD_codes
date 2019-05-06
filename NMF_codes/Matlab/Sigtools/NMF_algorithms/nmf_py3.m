function [W, H, reconstructError] = nmf_py3(X, params, varargin)


numberOfSources = params.numberOfSources;
init = params.init;
betaDivergence = params.betaDivergence;
numberOfIterations = params.numberOfIterations;
tolChange = params.tolChange;
tolError = params.tolError;

sqrteps = sqrt(eps);

betaDivergenceAux = beta_loss_to_float(betaDivergence);
[n_samples, n_features] = size(X);

switch init
    case 'random'
        avg = sqrt(mean(X(:)) / numberOfSources);
        H0 = avg * randn(numberOfSources, n_features); 
        W0 = avg * randn(n_samples, numberOfSources);
        H0 = abs(H0);
        W0 = abs(W0);
    case 'custom'
        if nargin > 2 && strcmp(init, 'custom')
            W0 = varargin{1};
            H0 = varargin{2};
            flag = checkMatrices(W0, H0);
            if flag
                error('Input matrices have nonnegative elements');
            end
        else
            error('Input initial matrices');
        end
        if size(W0, 2) ~= numberOfSources || size(H0, 1) ~= numberOfSources
            error('Input matrices sizes are incorrect');
        end
        
    otherwise
        error('Initialization method not supported');
end

reconstructError = zeros(numberOfIterations, 1);

for i = 1:numberOfIterations
    
    if betaDivergenceAux == 2
        numW = X*H0';
        W = max(0, W0 .* (numW ./ (W0*(H0*H0') + eps(numW))));
        numH = W'*X;
        H = max(0, H0 .* (numH ./ ((W'*W)*H0 + eps(numH))));
        
        
    elseif betaDivergenceAux == 1
        numW = (X./(W0*H0 + eps(X)))*H0.';
        denW = (ones(size(W0, 1), size(H0, 2))*H0.') + eps(numW);
        
        W = W0 .* (numW ./ denW); 
        
        numH = W.'*(X./(W*H0 + eps(X)));
        denH = W.'*ones(size(W, 1), size(H0, 2)) + eps(numH);
        
        H = H0 .* (numH ./ denH);
    end
    
    reconstructError(i) = beta_divergence(X, W, H, betaDivergence, true);
    dw = max(max(abs(W - W0) / (sqrteps + max(max(abs(W0))))));
    dh = max(max(abs(H - H0) / (sqrteps + max(max(abs(H0))))));
    delta = max(dw, dh);
    
    if i > 1
        if delta <= tolChange
            break;
        elseif (reconstructError0 - reconstructError(i)) <= (tolError*max(1, reconstructError0))
            break;
        elseif i == numberOfIterations
            disp('Max iterations reached. Increase this number for better convergence.')
            break
        end
    end
    
    reconstructError0 = reconstructError(i);
    W0 = W;
    H0 = H;
end

