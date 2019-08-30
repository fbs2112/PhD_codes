function [output] = detection_eval(y, threshold, varargin)

output = false(length(y), 1);
output(y > threshold) = true;

if nargin > 3
    if strcmp(varargin{2}, 'invert')
        output = false(length(y), 1);
        output(y < threshold) = true;
    end
end

if nargin == 3
    if ~mod(varargin{1}, 2)
        error('Window length must be odd');
    end
    window_length = varargin{1};
    output = window_eval(output, window_length, @median).';
end
