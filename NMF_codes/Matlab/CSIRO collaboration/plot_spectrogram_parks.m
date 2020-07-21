clear;
clc;
close all;

addpath(['.' filesep 'data']);

load dataParks.mat;

fs = Fs;
nfft = 256;
nperseg = 256;
overlap = nperseg - 1;
hop_size = nperseg - overlap;
window = kaiser(nperseg, 5);
specType = 'power';

[Pxx, f, t] = spectrogram(parksSignal, window, overlap, nfft, fs, 'centered', specType);

figure;
surf(t*1e6, f/1e6, 10*log10(abs(Pxx)), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlabel('Time');
ylabel('Frequency [MHz]');
ax = gca;
set(ax, 'colormap', jet);

rmpath(['.' filesep 'data']);