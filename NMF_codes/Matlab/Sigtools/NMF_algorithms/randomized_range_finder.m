function [Q] = randomized_range_finder(A, size, n_iter, power_iteration_normalizer, random_state)

rng(random_state);

Q = randn(size(A,2), size);

if strcmp(power_iteration_normalizer, 'auto')
    if n_iter <= 2
        power_iteration_normalizer = 'none';
    else
        power_iteration_normalizer = 'LU';
    end
end

for i = 1:n_iter
    if strcmp(power_iteration_normalizer, 'none')
        Q = A*Q;
        Q = A.'* Q;
    elseif strcmp(power_iteration_normalizer, 'LU')
        [Q, ~, P] = lu(A*Q);
        Q = P*Q;
        [Q, P, ~]= lu(A.'*Q);
        Q = P*Q;
    elseif strcmp(power_iteration_normalizer, 'QR')
        [Q, ~] = qr(A*Q, 0);
        [Q, ~] = qr(A.'*Q, 0);
    end
end

    % Sample the range of A using by linear projection of Q
    % Extract an orthonormal basis
[Q, ~] = qr(A*Q, 0);