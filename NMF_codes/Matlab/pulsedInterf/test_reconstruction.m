clear;
clc;
close all;

fs = 32.768e6;

x = randn(1000, 1) + 1j*randn(1000, 1);

pxx = spectrogram(x, ones(64, 1), 64 - 1, 64, fs, 'centered', 'power');

testCOLA = iscola(ones(64, 1), 63);

if testCOLA
    xHat = istft(pxx,fs,'Window',ones(64, 1),'OverlapLength',63,'FFTLength',64);
    figure;
    plot(abs(x));
    figure;
    plot(abs(xHat), 'r');
    
else
    error('Reconstruction not permitted');
end


