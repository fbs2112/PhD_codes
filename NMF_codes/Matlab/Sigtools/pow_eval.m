function [power] = pow_eval(y)

y = y(:);
power = y'*y/length(y);