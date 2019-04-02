function [W, H, n_iter] = non_negative_factorization(X, W, H, n_components, init, beta_loss, tol,...
                               max_iter, alpha, l1_ratio, regularization, random_state)
                           

                           
if ~strcmp(init, 'custom')
    [W, H] = initialize_nmf(X, n_components, init, random_state);
end

[l1_reg_W, l1_reg_H, l2_reg_W, l2_reg_H] = compute_regularization(alpha, l1_ratio, regularization);

[W, H, n_iter] = fit_multiplicative_update(X, W, H, beta_loss, max_iter, tol, l1_reg_W, l1_reg_H,...
                                                  l2_reg_W, l2_reg_H);
                                              
                                              
if n_iter == max_iter && tol > 0
    warning("Maximum number of iteration reached. Increase it to improve convergence.");
end