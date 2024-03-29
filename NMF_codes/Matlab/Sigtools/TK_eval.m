function [output] = TK_eval(x)

a = zeros(size(x));

for i = 1:size(x,2)
    y = x(:,i);
    b = y.^2;
    aux = [y(1); y(1:end-1)];
    aux2 = [y(2:end); y(end)];
    a(:,i) = b - aux2.*aux;
    a(1,i) = a(2,i);  %smoothing border effects
    a(end,i) = a(end-1,i); %smoothing border effects
end

output = a;