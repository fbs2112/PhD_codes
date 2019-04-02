function [W, H] = initialize_nmf(X, n_components, init, eps, random_state)

rng(random_state);

[n_samples, n_features] = size(X);

if strcmp(init, 'random')
    avg = sqrt(mean(X(:)) / n_components);
    H = avg * randn(n_components, n_features);
    W = avg * randn(n_samples, n_components);
    H = abs(H);
    W = abs(W);

elseif strcmp(init, 'nndsvd') || strcmp(init, 'nndsvda') || strcmp(init, 'nndsvdar')

    [U, S, V] = randomized_svd(X, n_components, random_state);
    W = zeros(size(U));
    H = zeros(size(V));    

    % The leading singular triplet is non-negative
    % so it can be used as is for initialization.
    W(:, 1) = sqrt(S(1,:))*abs(U(:, 1));
    H(1, :) = sqrt(S(1,:))*abs(V(1, :));

    for j = 2:n_components
        x = U(:, j);
        y = V(j, :);

        % extract positive and negative parts of column vectors
        x_p = max(x, 0);
        y_p = max(y, 0);
        x_n = abs(min(x, 0));
        y_n = abs(min(y, 0));

        % and their norms
        x_p_nrm = norm(x_p);
        y_p_nrm = norm(y_p);
        x_n_nrm = norm(x_n);
        y_n_nrm = norm(y_n);

        m_p = x_p_nrm * y_p_nrm;
        m_n = x_n_nrm * y_n_nrm;

        % choose update
        if m_p > m_n
            u = x_p / x_p_nrm;
            v = y_p / y_p_nrm;
            sigma = m_p;
        else
            u = x_n / x_n_nrm;
            v = y_n / y_n_nrm;
            sigma = m_n;
        end
        lbd = sqrt(S(j) * sigma);
        W(:, j) = lbd * u;
        H(j, :) = lbd * v;
    end
    
    W(W < eps) = 0;
    H(H < eps) = 0;
    
    if strcmp(init, 'nndsvda')
        avg = mean(X(:));
        W(W == 0) = avg;
        H(H == 0) = avg;
    end
    
    if strcmp(init, 'nndsvdar')
        avg = mean(X(:));
        W(W == 0) = abs(avg*randn(length(W(W == 0))) / 100);
        H(H == 0) = abs(avg*randn(length(H(H == 0))) / 100);
    end
        
else
    error('Initialization method not supported')
end


