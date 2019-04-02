function [output] = TK_filtering(x)

TK_result = TK_eval(x);

if isreal(TK_result)
    output = TK_result;
else
    output = TK_eval(real(x)) + TK_eval(imag(x));
end
   
   