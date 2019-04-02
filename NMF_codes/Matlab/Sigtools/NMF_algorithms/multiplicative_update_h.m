function [delta_H] = multiplicative_update_h(X, W, H, beta_loss, l1_reg_H, l2_reg_H)


if beta_loss == 2
    numerator = W.'*X;
    denominator = (W.'*W)* H;
    
elseif beta_loss == 1
    WH_safe_X = W*H;
    WH_safe_X_data = WH_safe_X;
    X_data = X;
    WH_safe_X_data = X_data./WH_safe_X_data;
    WH_safe_X = WH_safe_X_data;
    numerator = W.'*WH_safe_X;
    
    W_sum = sum(W, 1);  % shape(n_components, )
    W_sum(W_sum == 0) = 1;
    denominator = W_sum;
end

if l1_reg_H > 0
    denominator = denominator + l1_reg_H;
end
if l2_reg_H > 0
    denominator = denominator + l2_reg_H .* H;
    denominator(denominator == 0) = eps;
end

numerator = numerator./denominator;
delta_H = numerator;