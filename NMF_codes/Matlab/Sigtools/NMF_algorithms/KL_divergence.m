function [res] = KL_divergence(X, W, H)

       
WH = W*H;

div = X ./ WH;
res = (X).*log(div);
res = res + WH - X;
res = sum(res(:));
