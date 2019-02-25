clear;
clc;
close all;

addpath(['.' filesep 'Sigtools' filesep])

figPath = './figs-02-14/';
fs = 32.768e6;
nfft = 19;
nperseg = 19;
beta_loss = 'kullback-leibler';
save_fig = 0;

exp_name = 'TK';    
f0 = 0;
secondsOfData = 8.62e-6;
numberOfSamples = secondsOfData*fs;
totalSamples = 4096;
bandwidth = 1e6;
JNRVector = [-10, -5, 0, 10];
JNRVector = 10;


t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;

signal1 = exp(1j*2*pi*f.*t).';    
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
powSignal = pow_eval(signal1);
signal1 = [zeros(totalSamples - length(signal1), 1); signal1];
signalLength = length(signal1);
noise = randn(signalLength, 1) + 1j*randn(signalLength, 1);
powNoise = pow_eval(noise);

for i = 1:length(JNRVector)
    
    powAux = powSignal/JNRVector(i);
    noise2 = noise*sqrt(powAux/powNoise);
    data = signal1 + noise2;
    [PxxAux, f, t] = spectrogram(data, hann(nfft), nfft-1, nfft, fs, 'centered', 'power');
    inputNMF = abs(PxxAux).^2;
    specdB = 10*log10(inputNMF);
%     spectrogram(data, hann(nfft), nfft-1, nfft, fs, 'centered', 'power');

%     figure;
%     imshow(PxxAux)
end

rmpath(['.' filesep 'Sigtools' filesep])