function [U, s, V] = randomized_svd(M, n_components, n_oversamples, n_iter, power_iteration_normalizer, transpose, flip_sign, random_state)

n_random = n_components + n_oversamples;

[n_samples, n_features] = size(M);

if strcmp(n_iter, 'auto')
    if n_components < .1 * min(size(M))
        n_iter = 7;
    else
        n_iter = 4;
    end
end

if strcmp(transpose, 'auto')
    transpose = n_samples < n_features;
end

if transpose
    M = M.';
end

Q = randomized_range_finder(M, n_random, n_iter, power_iteration_normalizer, random_state);
B = Q.'*T;

[Uhat, s, V] = svd(B, 'econ');
U = Q*Uhat;

s = diag(s);

if flip_sign
    if ~transpose
       [U, V] = svd_flip(U, V, true(1));
    else
        [U, V] = svd_flip(U, V, false(1)); %see this function
    end
end

if transpose
    V = V(1:n_components+1, :).';
    s = s(1:n_components+1);
    U = U(:,1:n_components+1).';
else
    aux = V;
    V = U(:,1:n_components);
    U = aux(1:n_components,:);
end
    