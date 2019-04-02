function [res] = beta_divergence(X, W, H, beta, square_root)

beta = beta_loss_to_float(beta);

if beta == 2
    res = squared_norm(X - W*H) / 2;
end
       
WH = W*H;
WHaux = WH.';
WH_data = WHaux(:);
Xaux = X.';
X_data = Xaux(:);
indices = X_data > eps;
WH_data = WH_data(indices);
X_data = X_data(indices);
        
WH_data(WH_data == 0) = eps;

% generalized Kullback-Leibler divergence
if beta == 1
    sum_WH = sum(W, 1)*sum(H, 2);
    div = X_data ./ WH_data;
    res = (X_data).'*log(div);
    res = res + sum_WH - sum(X_data);
end

if square_root
    res = sqrt(res * 2);
end