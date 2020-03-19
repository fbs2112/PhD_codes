function [quantSignal, mse] = quantise_gps(signal, nbits, varargin)

Ag = 1;

if nargin == 3
    noiseVar = varargin{1};
    if isreal(signal)
        noiseStd = sqrt(noiseVar);
    else
        noiseStd = sqrt(noiseVar / 2);
    end
    switch nbits
        case 2
            Ag = 1 / noiseStd;
        case 3
            Ag = 1.71 / noiseStd;
        case 4
            Ag = 2.98/noiseStd;
        case 5
            Ag = 5.315/noiseStd;
        otherwise
            Ag = 1;
            fprintf('Warning, AGC gain not defined for nbits > 5\n');
            fprintf('AGC gain set to 1');
    end
end

signal = signal * Ag;
maxValue = (2^nbits - 1);
minValue = -maxValue;

q = (maxValue-minValue)/(2^nbits-1);
realPart = min(max(real(signal), minValue), maxValue);
realPart = minValue + round((realPart-minValue)/q)*q;
imagPart = 0;

if ~isreal(signal)
    imagPart = min(max(imag(signal), minValue), maxValue);
    imagPart = minValue + round((imagPart-minValue)/q)*q;
end

quantSignal = realPart + 1j*imagPart;
mse = norm(signal - quantSignal).^2;