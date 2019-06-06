
function [window] = window_eval(y, window_length, function_eval, varargin)

y = y(:);
hop_size = 1;

if nargin > 3
    hop_size = varargin{1};
end

aux = (length(y) - window_length)/hop_size;

if aux < 0
    error('Window size and hop size are incompatible.')
end

inputLength = length(y);
y = [y(1)*ones(window_length-1, 1); y ; y(end)*ones(window_length-1, 1)]; % padding with the first and last elements of y
yAux = buffer(y, window_length, window_length - hop_size, 'nodelay');
window = function_eval(yAux(:, 1:inputLength));
