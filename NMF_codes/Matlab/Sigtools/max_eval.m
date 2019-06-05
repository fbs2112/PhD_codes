function [output] = max_eval(x, varargin)

if nargin == 1
    axis = 1;
else
    axis = varargin{1};
end

output = max(x, [], axis); 