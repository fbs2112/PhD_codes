function [flag] = checkMatrices(W, H)

W = W(:);
H = H(:);

flag = false(1);

if any(W < 0) || any(H < 0)
    flag = true(1);
end
