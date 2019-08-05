function [power] = pow_eval(y)

power = mean(abs(y).^2, 1);
