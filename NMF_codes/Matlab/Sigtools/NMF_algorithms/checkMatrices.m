%This function checks whether the input matrices W and H have negative
%elements. If that is the case, it returns true, and false otherwise
%
%This code was created by Felipe Barboza da Silva
% Copyright (c) School of Engineering. Macquarie University - Australia.
% All rights reserved. 2020.
% DISCLAIMER:
%     This code is strictly private, confidential and personal to its recipients
%     and should not be copied, distributed or reproduced in whole or in part,
%     nor passed to any third party.
% ** Do not duplicate or distribute without written permission from the owners. **

function [flag] = checkMatrices(W, H)

W = W(:);
H = H(:);

flag = false(1);

if any(W < 0) || any(H < 0)
    flag = true(1);
end
