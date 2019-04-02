function [u, v, signs] = svd_flip(u, v, u_based_decision)


if u_based_decision
    % columns of u, rows of v
    [~, max_abs_cols] = max(abs(u), [], 1);
%     for i = 1:
    signs = sign(u(max_abs_cols, 1:size(u, 2))); % this line in python brings different results
    u = u.*signs;
    v = v.*signs;
else
    % rows of v, columns of u
    [~, max_abs_rows] = np.argmax(abs(v), [], 2);
    signs = sign(v(1:size(v, 1), max_abs_rows));
    u = u.*signs;
    v = v.*signs;
end