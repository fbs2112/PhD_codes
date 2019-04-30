function [output] = detection_eval(y, threshold, varargin)

output = false(length(y), 1);
output(y > threshold) = true;

if nargin > 2
    if ~mod(varargin{1}, 2)
        error('Window length must be odd');
    end
    window_length = varargin{1};
    output = window_eval(output, window_length, @median).';
end
