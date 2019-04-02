function [beta_loss] = beta_loss_to_float(beta_loss)

switch beta_loss
    
    case 'frobenius'
        beta_loss = 2;
    case 'kullback-leibler'
        beta_loss = 1;
        
    otherwise
        error('Beta loss not supported');
end
        