function [output] = similarity_eval(x, y, similarity, varargin)

x = x(:);
y = y(:);

if nargin > 3
    if varargin{2}
         x = x - mean(x);
         x =  x*sqrt(1/var(x));
         y = y - mean(y);
         y = y*sqrt(1/var(y));
    end
end

switch similarity
    case 'inner'
        output = (x.'*y);
        if strcmp(varargin{1}, 'normalized')
            output = output/(norm(x)*norm(y) + eps);
        end
        
    case 'gaussian'
        if nargin < 3
            error('Please choose a standard deviation value');
        end
        
        sigma = varargin{1};
        output =  exp(-(norm(x - y)^2)/(2*sigma^2));
        
    case 'euc'       
        output =  1/(norm(x - y) + eps);
        
    case 'corr'
        output =  sum(xcorr(x, y));
        
        
    case 'polinomial'
        if nargin < 3
            error('Please choose the polinomial degree');
        end
        disp('Still needs to be implemented');
        output = [];
        
    otherwise
        error('Please type inner, gausssian or polinomial as similarity function');
        
end