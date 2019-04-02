
function [window] = window_eval(y, window_length, function_eval, varargin)

hop_size = 1;

y(y==0) = eps;

if nargin > 3
    hop_size = varargin{1};
end

aux = (length(y) - window_length)/hop_size;

 if aux < 0
    error('Window size and hop size are incompatible.')
end

yAux = buffer(y, window_length, window_length - hop_size, 'nodelay');

lastColumn = yAux(:, end);

idxNonzero = find(lastColumn, 1, 'last');
idxZero = find(~lastColumn);
yAux(idxZero,end) = lastColumn(idxNonzero);
% number_of_windows = ceil(aux + 1);
window = function_eval(yAux, 1);
