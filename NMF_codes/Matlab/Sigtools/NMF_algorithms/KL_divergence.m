%This function evaluates the Kullback-leibler divergence given the matrices
%X, W and H.
%
%This code was created by Felipe Barboza da Silva
% Copyright (c) School of Engineering. Macquarie University - Australia.
% All rights reserved. 2020.
% DISCLAIMER:
%     This code is strictly private, confidential and personal to its recipients
%     and should not be copied, distributed or reproduced in whole or in part,
%     nor passed to any third party.
% ** Do not duplicate or distribute without written permission from the owners. **

function [res] = KL_divergence(X, W, H, varargin)
       
WH = W*H;

div = X ./ WH;
res = (X).*log(div);
res = res + WH - X;
res = sum(res(:));

if nargin == 4
    mu = varargin{1};
    res = res + mu(1)*norm(W, 1) + mu(2)*norm(H, 1);
end