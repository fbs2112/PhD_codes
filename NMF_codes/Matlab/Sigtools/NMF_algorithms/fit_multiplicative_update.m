function [W, H, index] = fit_multiplicative_update(X, W, H, beta_loss, max_iter, tol, l1_reg_W, l1_reg_H, l2_reg_W, l2_reg_H)

beta_loss = beta_loss_to_float(beta_loss);

error_at_init = beta_divergence(X, W, H, beta_loss, true(1));

previous_error = error_at_init;

H_sum = [];
HHt = [];
XHt = [];

for index = 1:max_iter
    [delta_W, H_sum, HHt, XHt] = multiplicative_update_w(X, W, H, beta_loss, l1_reg_W, l2_reg_W, H_sum, HHt, XHt);
    W = W.*delta_W;
    
    delta_H = multiplicative_update_h(X, W, H, beta_loss, l1_reg_H, l2_reg_H);
    H = H.*delta_H;
    H_sum = [];
    HHt = [];
    XHt = [];
    
    if tol > 0 && mod(index, 10) == 0
        error = beta_divergence(X, W, H, beta_loss, true(1));
            
        if ((previous_error - error) / error_at_init) < tol
            break
            previous_error = error;
        end
    end
end