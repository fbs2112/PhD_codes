clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['.' filesep 'data']);

nfft = 256;
nperseg = 256;
overlap = nperseg/2;
% window = kaiser(nperseg, 5);
window = ones(nperseg, 1);
specType = 'power';

freqVector = linspace(-64, -56, nfft)*1e6;

for i = 2:198
    load(['.' filesep 'data' filesep 'dataParkes_' num2str(i) '.mat']);
    fs = Fs;
    [Pxx, f, t] = spectrogram(parkesSignal, window, overlap, nfft, fs, specType);
    figure;
    surf(t*1e3, f/1e6, 10*log10(abs(Pxx)), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlim([min(t) max(t)]*1e3);
%     ylim([min(f)/1e6 -56]);
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    ax = gca;
    set(ax, 'colormap', jet);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    prompt = 'Type 1 for interference present or 0 otherwise\n';
    interferenceDetFlag(i) = input(prompt);
    close all;
end
rmpath(['.' filesep 'data']);