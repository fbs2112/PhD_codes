function [delta_W, H_sum, HHt, XHt] = multiplicative_update_w(X, W, H, beta_loss, l1_reg_W, l2_reg_W, H_sum, HHt, XHt)

if beta_loss == 2
    XHt = X*H.';
    numerator = XHt;
    HHt = H*H.';
    denominator = W*HHt;
    H_sum = [];
elseif beta_loss == 1
    HHt = [];
    XHt = [];
    WH_safe_X = W*H;
    WH_safe_X_data = WH_safe_X;
    X_data = X;
    WH_safe_X_data = X_data./WH_safe_X_data;
    WH_safe_X = WH_safe_X_data; %these attributions come from Python
    
    numerator = WH_safe_X*H.';
    H_sum = sum(H, 2);
    denominator = H_sum.';
else
    error('Divergence not supported');
end

% Add L1 and L2 regularization
if l1_reg_W > 0
    denominator = denominator + l1_reg_W;
end
if l2_reg_W > 0
    denominator = denominator + l2_reg_W.*W;
end
denominator(denominator == 0) = eps;

numerator = numerator./denominator;
delta_W = numerator;


    
    
