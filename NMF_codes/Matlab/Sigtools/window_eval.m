
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

y = [y(1)*ones(window_length-1, 1); y]; % padding with the first element of y
y(y==0) = eps;

yAux = buffer(y, window_length, window_length - hop_size, 'nodelay');

lastColumn = yAux(:, end);

idxNonzero = find(lastColumn, 1, 'last');
idxZero = find(~lastColumn);
yAux(idxZero,end) = lastColumn(idxNonzero);
window = function_eval(yAux, 1);
window(window<=eps) = 0;
