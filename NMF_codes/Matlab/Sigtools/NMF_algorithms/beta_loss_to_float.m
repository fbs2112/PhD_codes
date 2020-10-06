%This function assigns an integer variable according to the input string
%beta_loss. If beta_loss is 'euclidean', it returns 2. If beta_loss is
%'kullback-leibler, it returns 1.
%
%This code was created by Felipe Barboza da Silva
% Copyright (c) School of Engineering. Macquarie University - Australia.
% All rights reserved. 2020.
% DISCLAIMER:
%     This code is strictly private, confidential and personal to its recipients
%     and should not be copied, distributed or reproduced in whole or in part,
%     nor passed to any third party.
% ** Do not duplicate or distribute without written permission from the owners. **

function [beta_loss] = beta_loss_to_float(beta_loss)

switch beta_loss
    
    case 'euclidean'
        beta_loss = 2;
    case 'kullback-leibler'
        beta_loss = 1;
        
    otherwise
        error('Beta loss not supported');
end
        