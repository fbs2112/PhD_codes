function [res] = beta_divergence_v2(X, W, H, beta, varargin)

sigma = W*H;

if beta == 2
    res = norm((X - sigma) / 2, 'fro');
elseif beta == 1
    res = X.*log(X./sigma) + sigma - X;
    res = sum(res(:));
elseif beta == 0
    res = X./sigma - log(X./sigma) - 1;
    res = sum(res(:));
end
       
if nargin == 4
    mu = varargin{1};
    res = res + mu(1)*norm(W(:), 1) + mu(2)*norm(H(:), 1);
end