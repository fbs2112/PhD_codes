function [quantSignal, mse] = quantise_gps(signal, nbits)

minValue = min(signal);
maxValue = max(signal);

imagPart = 0;
q = (maxValue-minValue)/(2^nbits-1);
realPart = min(max(real(signal), minValue), maxValue);
realPart = minValue + round((realPart-minValue)/q)*q;

if ~isreal(signal)
    imagPart = min(max(imag(signal), minValue), maxValue);
    imagPart = minValue + round((imagPart-minValue)/q)*q;
end

quantSignal = realPart + 1j*imagPart;
mse = norm(signal - quantSignal).^2;